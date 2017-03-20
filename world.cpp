#include <random>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <chrono>

#include "world.h"
#include "patch.h"
#include "individual.h"

World::World(int idWorld, int NPatch, double delta, double c, int typeMut, double mu, double sigmaZ, int Kmin, int Kmax, int sigmaK,
             double Pmin, double Pmax, double sigmaP, double sInit, double dInit, int NGen, int genReport)
{
    int i = 0;

    this->NPatch = NPatch;

    this->delta = delta;
    this->c = c;

    this->typeMut = (distrMut)typeMut;
    this->mu = mu;
    this->sigmaZ = sigmaZ;

    this->NGen = NGen;

    report.open ("report_" + std::to_string(idWorld) + ".txt");
    this->genReport = genReport;

    /* Les lignes suivantes permettent de réserver
    de la mémoire pour éviter les réallocations
    qui peuvent diminuer les performances. */
    patches.reserve(NPatch);
    juveniles[0].reserve(Kmax);
    juveniles[1].reserve(Kmax);

    for(i=0; i<NPatch; i++)
    {
        patches.emplace_back(distr(Pmin, Pmax, sigmaP, i), distr(Kmin, Kmax, sigmaK, i), sInit, dInit);
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed (seed);

    writeHeader(Kmin, Kmax, sigmaK, Pmin, Pmax, sigmaP);
}

double World::distr(double minVal, double maxVal, double sigma, int posPatch)
{
    return (minVal + ((maxVal - minVal) * exp( - ((posPatch-(NPatch/2))*(posPatch-(NPatch/2))) / (2*sigma))));
}

void World::run(int idWorld)
{
    int i = 0, j = 0, progress = 0;

    std::cout << "Progression du monde " << idWorld << " :" << std::endl;
    printProgress(0);

    for(i=0; i<=NGen; i++)
    {
        if(i%genReport == 0) {writeReport(i);}

        for(j=0; j<NPatch; j++) {createNextGen(j);}

        for(j=0; j<NPatch; j++) {patches[j].isPollenized();}

        /* Indique la progression à l'écran */
        if (i*78/NGen > progress)
        {
            progress = i*78/NGen;
            printProgress(progress);
        }
    }
}

void World::createNextGen (int idPatch)
{
    int i = 0;

    /* Vecteur qui contient toutes les mères possibles.
       Elles sont numérotées de 0 à n. */
    std::vector<int> mother;

    /* Vecteur qui contient toutes les pressions
       pour un patch (dispersantes des voisins
       et résidentes du patch local).
       Les valeurs paires sont issues d'autof.
       Les valeurs impaires sont issues d'allof. */
    std::vector<double> press;

    /* Vecteur qui permet de savoir d'où vient la mère retenue */
    std::vector<int> fromPatch;

    /* Les lignes suivantes permettent de réserver
    de la mémoire pour éviter les réallocations
    qui peuvent diminuer les performances. */
    mother.reserve(6*patches[idPatch].K +1);
    press.reserve(6*patches[idPatch].K);
    fromPatch.reserve(6*patches[idPatch].K);

    if (idPatch != 0)
    {
        patches[idPatch - 1].getPression(delta, c, true, press);

        /* On peut vider le vecteur car il n'est plus utile. */
        clear_and_freeVector(patches[idPatch - 1].dispSeeds);

        fromPatch.insert(fromPatch.end(), 2*patches[idPatch - 1].K, idPatch - 1);
    }

    patches[idPatch].getPression(delta, c, false, press);

    fromPatch.insert(fromPatch.end(), 2*patches[idPatch].K, idPatch);

    if (idPatch != NPatch - 1)
    {
        patches[idPatch + 1].getPression(delta, c, true, press);

        fromPatch.insert(fromPatch.end(), 2*patches[idPatch + 1].K, idPatch + 1);
    }

    for (i=0; i<(int)press.size() + 1; i++)
    {
        mother.push_back(i);
    }

    /* Objet qui permet de générer des nombres aléatoires pondérés. */
    std::piecewise_constant_distribution<double> weighted (mother.begin(), mother.end(), press.begin());

    for(i=0; i<patches[idPatch].K; i++)
    {
        int chosenMother = weighted(generator);
        int patchMother = fromPatch[chosenMother];

        bool autof = false;
        if (chosenMother%2 == 0) {autof = true;}

        if (idPatch == 0)
        {
            if (fromPatch[chosenMother] == idPatch + 1)
            {
                chosenMother = chosenMother - 2*patches[idPatch].K;
            }
        }

        else
        {
            if (fromPatch[chosenMother] == idPatch)
            {
                chosenMother = chosenMother - 2*patches[idPatch - 1].K;
            }

            if (fromPatch[chosenMother] == idPatch + 1)
            {
                chosenMother = chosenMother - 2*patches[idPatch - 1].K - 2*patches[idPatch].K;
            }
        }

        chosenMother = chosenMother/2;

        /* Pour les patchs pairs, on met la nouvelle génération dans le 1er vecteur.
        Pour les patchs impairs, dans le 2nd. */
        newInd(idPatch%2, patchMother, chosenMother, autof);
        mutation(juveniles[idPatch%2][i]);
    }

    /* Au premier patch, rien à faire. */
    if(idPatch != 0)
    {
        patches[idPatch - 1].population = juveniles[(idPatch-1)%2];
        juveniles[(idPatch-1)%2].clear();
    }

    /* Au dernier patch, on remplace la génération. */
    if(idPatch == NPatch - 1)
    {
        patches[idPatch].population = juveniles[idPatch%2];
        juveniles[idPatch%2].clear();
    }
}

void World::newInd(int whr, int patchMother, int mother, bool autof)
{
    std::array<double,2> fatherTraits = {0,0};

    /* Issue d'autof */
    if(autof)
    {
        juveniles[whr].emplace_back(patches[patchMother].population[mother].s,
                              patches[patchMother].population[mother].d);
    }

    /* Sinon, on cherche un père */
    else
    {
        getFather(patchMother, mother, fatherTraits);
        juveniles[whr].emplace_back((patches[patchMother].population[mother].s + fatherTraits[0])/2,
                              (patches[patchMother].population[mother].d + fatherTraits[1])/2);
    }
}

void World::mutation(Individual& IndToMutate)
{
    /* Pour «retenir» quel trait doit muter
       false: d     true: s */
    bool sWasChosen = false;
    double t = IndToMutate.d;

    std::uniform_real_distribution<double> unif(0, 1);

    /* Y a-t-il mutation ? */
    if(unif(generator) < mu)
    {
        /* On choisit quel trait mute */
        if(unif(generator) >= 0.5)
        {
            t = IndToMutate.s;
            sWasChosen = true;
        }

        switch(typeMut)
        {
            case gaussian:
                t = gaussMutation(t);
                break;

            case uniform:
                t = unifMutation(t);
                break;
        }

        if (sWasChosen) {IndToMutate.s = t;}
        else {IndToMutate.d = t;}
    }
}

double World::gaussMutation (double t)
{
    std::normal_distribution<double> gauss(0,sigmaZ);
    double deltaMu = gauss(generator);

    return t*exp(deltaMu)/((exp(deltaMu) - 1)*t + 1);
}

double World::unifMutation (double t)
{
    double lowerBound = t - sigmaZ;
    double upperBound = t + sigmaZ;

    if(lowerBound < 0) {lowerBound = 0;}
    if(upperBound > 1) {upperBound = 1;}

    std::uniform_real_distribution<double> unif(lowerBound, upperBound);

    return unif(generator);
}

void World::getFather(int patchMother, int mother, std::array<double,2>& fatherTraits)
{
    int father = 0;

    std::uniform_int_distribution<int> unif(0, patches[patchMother].K-1);

    do {father = unif(generator);}
    while (father == mother); //Pas de pseudo allofécondation

    fatherTraits[0] = patches[patchMother].population[father].s;
    fatherTraits[1] = patches[patchMother].population[father].d;
}

void World::clear_and_freeVector(std::vector<double>& toClear)
{
    toClear.clear();
    std::vector<double>().swap(toClear);
}

void World::printProgress(int progress)
{
    int i = 0;
    std::cout << "[";
    for (i=0; i<78; i++)
    {
        if(i<progress) {std::cout << "#";}
        else {std::cout << ".";}
    }
    std::cout << "]" << std::endl;
}

void World::writeHeader(int Kmin, int Kmax, int sigmaK, double Pmin, double Pmax, double sigmaP)
{
    report << "Nombre de patchs=" << NPatch << std::endl;
    report << "Delta=" << delta << " c=" << c << std::endl;
    report << "Loi pour la mutation:" << typeMut << " mu=" << mu << " sigmaZ=" << sigmaZ << std::endl;
    report << "Kmin=" << Kmin << " Kmax=" << Kmax << " SigmaK=" << sigmaK << std::endl;
    report << "Pmin=" << Pmin << " Pmax=" << Pmax << " SigmaP=" << sigmaP << std::endl;
    report << "Gen\tPatch\tInd\ts\td" << std::endl;
}

void World::writeReport(int gen)
{
    int i = 0, j = 0;


    for(j=0; j<NPatch; j++)
    {
        for(i=0; i<patches[j].K; i++)
        {
            report << gen << '\t';
            report << j << '\t';
            report << i << '\t';
            report << patches[j].population[i].s << '\t';
            report << patches[j].population[i].d << std::endl;
        }
    }


    if (gen == NGen) {report.close();}
}

