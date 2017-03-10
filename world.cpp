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

World::World(int id, int NPatch, double delta, double c, int typeMut, double mu, double sigmaZ, int Kmin, int Kmax, int sigmaK,
             double Pmin, double Pmax, double sigmaP, double sInit, double dInit, int NGen, int genReport)
{
    int i = 0;

    this->id = id;

    this->NPatch = NPatch;

    this->delta = delta;
    this->c = c;

    this->typeMut = (distrMut)typeMut;
    this->mu = mu;
    this->sigmaZ = sigmaZ;

    this->Kmin = Kmin;
    this->Kmax = Kmax;
    this->sigmaK = sigmaK;

    this->Pmin = Pmin;
    this->Pmax = Pmax;
    this->sigmaP = sigmaP;

    this->NGen = NGen;

    report.open ("report_" + std::to_string(id) + ".txt");
    this->genReport = genReport;

    for(i=0; i<NPatch; i++)
    {
        patches.emplace_back(distr(Pmin, Pmax, sigmaP, i), distr(Kmin, Kmax, sigmaK, i), sInit, dInit);
    }

    writeHeader();
}

double World::distr(double minVal, double maxVal, double sigma, int i)
{
    return (minVal + ((maxVal - minVal) * exp( - ((i-(NPatch/2))*(i-(NPatch/2))) / (2*sigma))));
}

void World::run(void)
{
    int i = 0, j = 0, progress = 0;

    std::cout << "Progression du monde " << id << " :" << std::endl;
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

void World::createNextGen (int id)
{
    int i = 0;

    /* Seed basée sur l'horloge pour le générateur de nombres aléatoires */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    /* Générateur de nombres pseudo-aléatoires */
    std::default_random_engine generator (seed);

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

    if (id != 0)
    {
        patches[id - 1].getPression(delta, c, true, press);

        for (i=0; i<2*patches[id - 1].K; i++)
        {
            fromPatch.push_back(id - 1);
        }
    }

    patches[id].getPression(delta, c, false, press);

    for (i=0; i<2*patches[id].K; i++)
    {
        fromPatch.push_back(id);
    }

    if (id != NPatch - 1)
    {
        patches[id + 1].getPression(delta, c, true, press);

        for (i=0; i<2*patches[id + 1].K; i++)
        {
            fromPatch.push_back(id + 1);
        }
    }

    for (i=0; i<(int)press.size() + 1; i++)
    {
        mother.push_back(i);
    }

    /* Objet qui permet de générer des nombres aléatoires pondérés. */
    std::piecewise_constant_distribution<double> weighted (mother.begin(), mother.end(), press.begin());

    for(i=0; i<patches[id].K; i++)
    {
        int chosenMother = weighted(generator);
        int patchMother = fromPatch[chosenMother];

        bool autof = false;
        if (chosenMother%2 == 0) {autof = true;}

        if (id == 0)
        {
            if (fromPatch[chosenMother] == id + 1)
            {
                chosenMother = chosenMother - 2*patches[id].K;
            }
        }

        else
        {
            if (fromPatch[chosenMother] == id)
            {
                chosenMother = chosenMother - 2*patches[id - 1].K;
            }

            if (fromPatch[chosenMother] == id + 1)
            {
                chosenMother = chosenMother - 2*patches[id - 1].K - 2*patches[id].K;
            }
        }

        chosenMother = chosenMother/2;

        newInd(id%2, patchMother, chosenMother, autof);
        tmp[id%2][i].mutation(mu, sigmaZ, typeMut);
    }

    /* Au premier patch, rien à faire */
    if(id != 0)
    {
        patches[id - 1].population = tmp[(id-1)%2];
        tmp[(id-1)%2].clear();
    }

    /* Au dernier patch, on remplace la génération */
    if(id == NPatch - 1)
    {
        patches[id].population = tmp[(id)%2];
        tmp[(id)%2].clear();
    }
}

void World::newInd(int whr, int patchMother, int mother, bool autof)
{
    std::array<double,2> fatherTraits = {0,0};

    /* Issue d'autof */
    if(autof)
    {
        tmp[whr].emplace_back(patches[patchMother].population[mother].s,
                              patches[patchMother].population[mother].d);
    }

    /* Sinon, on cherche un père */
    else
    {
        getFather(patchMother, mother, fatherTraits);
        tmp[whr].emplace_back((patches[patchMother].population[mother].s + fatherTraits[0])/2,
                              (patches[patchMother].population[mother].d + fatherTraits[1])/2);
    }
}

void World::getFather(int idPatch, int mother, std::array<double,2>& fatherTraits)
{
    int father = 0;

    /* Seed basée sur l'horloge pour le générateur de nombres aléatoires */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    /* Générateur de nombres pseudo-aléatoires */
    std::default_random_engine generator (seed);

    std::uniform_int_distribution<int> unif(0, patches[idPatch].K-1);

    do {father = unif(generator);}
    while (father == mother); //Pas de pseudo allofécondation

    fatherTraits[0] = patches[idPatch].population[father].s;
    fatherTraits[1] = patches[idPatch].population[father].d;
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

void World::writeHeader(void)
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

