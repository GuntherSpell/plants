#include <random>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <chrono>
#include <algorithm>

#include "world.h"
#include "patch.h"
#include "individual.h"

World::World(int idWorld, int NPatch, double delta, double c, bool relationshipIsManaged,
             int typeMut, double mu, double sigmaZ, double d_s_relativeMutation, int Kmin, int Kmax, int sigmaK,
             double Pmin, double Pmax, double sigmaP, double sInit, double dInit, int NGen, int genReport)
{
    int i = 0, j = 0;

    this->NPatch = NPatch;

    this->delta = delta;
    this->c = c;

    this->relationshipIsManaged = relationshipIsManaged;

    this->typeMut = distrMut(typeMut);
    this->mu = mu;
    this->sigmaZ = sigmaZ;
    this->d_s_relativeMutation = d_s_relativeMutation;

    this->NGen = NGen;
    genCount = 0;

    report.open ("report_" + std::to_string(idWorld) + ".txt");
    this->genReport = genReport;

    /* Les lignes suivantes permettent de réserver
    de la mémoire pour éviter les réallocations
    qui peuvent diminuer les performances. */
    patches.reserve(NPatch);
    juveniles[0].reserve(Kmax);
    juveniles[1].reserve(Kmax);

    int Ktot = 0;

    for(i=0; i<NPatch; i++)
    {
        patches.emplace_back(distr(Pmin, Pmax, sigmaP, i), distr(Kmin, Kmax, sigmaK, i), sInit, dInit, Ktot);
        Ktot += patches[i].K;
    }

    globalPop.reserve(Ktot);
    for(i=0; i<NPatch; i++)
    {
        for(j=0; j<patches[i].K; j++)
        {
            IndividualPosition InfoToAdd;   //
            InfoToAdd.patch = i;            // On ne peut pas ajouter l'info dans le vecteur sans la construire avant.
            InfoToAdd.posInPatch = j;       //

            globalPop.push_back(InfoToAdd);
        }
    }

    if(relationshipIsManaged)
    {
        relationship[0].reserve(Ktot);
        relationship[1].reserve(Ktot);
        fathers.reserve(Ktot);
        mothers.reserve(Ktot);

        /* On part d'individus non apparentés. */
        for(i=0; i<Ktot; i++)
        {
            /* Pour chacune des deux matrices, on construit Ktot lignes remplies de i+1 zéros. */
            relationship[0].emplace_back(i + 1);
            relationship[1].emplace_back(i + 1);
        }
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed (seed);

    writeHeader(Kmin, Kmax, sigmaK, Pmin, Pmax, sigmaP);
    writeReport();
}

double World::distr(double minVal, double maxVal, double sigma, int posPatch)
{
    return (minVal + ((maxVal - minVal) * exp( - ((posPatch-(NPatch/2))*(posPatch-(NPatch/2))) / (2*sigma*sigma))));
}

void World::run(int idWorld)
{
    int i = 0, progress = 0, reportCount = 0;

    std::cout << "Progression du monde " << idWorld << " :" << std::endl;
    printProgress(0);

    for(genCount=0; genCount<=NGen; genCount++)
    {
        for(i=0; i<NPatch; i++)
        {
            patches[i].pollenized = redefinePollination(patches[i].p);
        }

        for(i=0; i<NPatch; i++) {createNextGen(i);}

        if(genCount%genReport == 0)
        {
            writeReport();

            if(reportCount < 2)
            {
                for(i=0; i<NPatch; i++)
                {
                    int j = 0;

                    std::array<double, 2> mean_var = {0,0};

                    std::vector<double> list_of_d;
                    std::vector<double> list_of_s;

                    for(j=0; j<patches[i].K; j++)
                    {
                        list_of_d.push_back(patches[i].population[j].d);
                        list_of_s.push_back(patches[i].population[j].s);
                    }

                    patches[i].calc_mean_var(list_of_d, mean_var);
                    patches[i].previous_d_means[reportCount%2] = mean_var[0];
                    patches[i].previous_d_vars[reportCount%2] = mean_var[1];

                    patches[i].calc_mean_var(list_of_s, mean_var);
                    patches[i].previous_s_means[reportCount%2] = mean_var[0];
                    patches[i].previous_s_vars[reportCount%2] = mean_var[1];
                }
            }

            else
            {
                for(i=0; i<NPatch; i++)
                {
                    int j = 0;
                    std::array<double, 2> mean_var = {0,0};

                    if(!patches[i].d_hasConverged)
                    {
                        std::vector<double> list_of_d;
                        for(j=0; j<patches[i].K; j++)
                        {
                            list_of_d.push_back(patches[i].population[j].d);
                        }

                        patches[i].calc_mean_var(list_of_d, mean_var);

                        patches[i].d_hasConverged = patches[i].check_stats(patches[i].previous_d_means, patches[i].previous_d_vars, mean_var);

                        patches[i].previous_d_means[reportCount%2] = mean_var[0];
                        patches[i].previous_d_vars[reportCount%2] = mean_var[1];
                    }

                    if(!patches[i].s_hasConverged)
                    {
                        std::vector<double> list_of_s;
                        for(j=0; j<patches[i].K; j++)
                        {
                            list_of_s.push_back(patches[i].population[j].s);
                        }

                        patches[i].calc_mean_var(list_of_s, mean_var);

                        patches[i].s_hasConverged = patches[i].check_stats(patches[i].previous_s_means, patches[i].previous_s_vars, mean_var);

                        patches[i].previous_s_means[reportCount%2] = mean_var[0];
                        patches[i].previous_s_vars[reportCount%2] = mean_var[1];
                    }
                }
            }

            reportCount ++;

            int sumOfConvergedPatches = 0;
            for(i=0; i<NPatch; i++)
            {
                if(patches[i].d_hasConverged && patches[i].s_hasConverged)
                {
                    sumOfConvergedPatches ++;
                }
            }

            if(sumOfConvergedPatches >NPatch*2/3) {return;}
        }

        if(relationshipIsManaged) {calcNewRelationships();}

        /* Indique la progression à l'écran */
        if (genCount*78/NGen > progress)
        {
            progress = genCount*78/NGen;
            printProgress(progress);
        }
    }
}

void World::createNextGen (int idPatch)
{
    int i = 0;

    /* Pour savoir la taille des vecteurs et ainsi réserver
    en avance la mémoire pour améliorer les performances. */
    int memoryToReserve = 2*patches[idPatch].K;

    /* Perment de savoir où commencer le tirage aléatoire des mères. */
    int firstMother = patches[idPatch].pos_of_first_ind;

    /* Vecteur qui contient toutes les mères possibles (allof et autof).
    Elles sont numérotées de firstMother à n.
    La 1ère mère fait de l'autof, la 2nde allof, 3ème autof, etc. */
    std::vector<int> mother;

    /* Vecteur qui contient toutes les pressions pour un patch
    (dispersantes des voisins et résidentes du patch local).
    Les valeurs paires sont issues d'autof.
    Les valeurs impaires sont issues d'allof. */
    std::vector<double> press;

    /* Selon la position du patch, il faut réserver plus ou moins de mémoire. */
    if(idPatch != 0)
    {
        memoryToReserve += 2*patches[idPatch - 1].K;
    }
    if(idPatch != NPatch - 1)
    {
        memoryToReserve += 2*patches[idPatch + 1].K;
    }

    mother.reserve(memoryToReserve + 1);
    press.reserve(memoryToReserve);

    if (idPatch != 0)
    {
        patches[idPatch - 1].getDispPress(delta, c, press);

        /* On peut vider le vecteur car il n'est plus utile. */
        clear_and_freeVector(patches[idPatch - 1].dispSeeds);

        /* Puisqu'on n'est pas tout à gauche, la première mère devient le premier individu du patch de gauche. */
        firstMother = patches[idPatch - 1].pos_of_first_ind;
    }

    patches[idPatch].getResidPress(delta, press);

    if (idPatch != NPatch - 1)
    {
        patches[idPatch + 1].getDispPress(delta, c, press);
    }

    for (i=firstMother; i<firstMother + int(press.size()) + 1; i++)
    {
        mother.push_back(i);
    }

    /* Objet qui permet de générer des nombres aléatoires pondérés. */
    std::piecewise_constant_distribution<double> weighted (mother.begin(), mother.end(), press.begin());

    /* Il faut savoir si les mères paires font de l'autof ou de l'allof. Cela dépend de firstMother. */

    for(i=0; i<patches[idPatch].K; i++)
    {
        int chosenMother = weighted(generator);

        /* Une mère fait de l'autof si elle a la même parité que la première mère. */
        bool autof = false;
        if (chosenMother%2 == firstMother%2)
        {
            autof = true;
        }

        chosenMother = (firstMother + chosenMother)/2;

        /* Pour les patchs pairs, on met la nouvelle génération dans le 1er vecteur.
        Pour les patchs impairs, dans le 2nd. */
        newInd(idPatch%2, chosenMother, autof);
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

void World::newInd(int whr, int mother, bool autof)
{
    double f = 0;

    /* On récupère la position relative de la mère dans son patch. */
    int patchMother = globalPop[mother].patch;
    int mother_PosInPatch = globalPop[mother].posInPatch;

    /* Issue d'autof */
    if(autof)
    {
        if(relationshipIsManaged)
        {
            f = 1/2 + patches[patchMother].population[mother_PosInPatch].f/2;
            mothers.push_back(mother);
            fathers.push_back(mother);
        }

        juveniles[whr].emplace_back(patches[patchMother].population[mother_PosInPatch].s,
                                    patches[patchMother].population[mother_PosInPatch].d, f);

    }

    /* Sinon, on cherche un père. */
    else
    {
        int father = getFather(patchMother, mother_PosInPatch);

        if(relationshipIsManaged)
        {
            mothers.push_back(mother);
            fathers.push_back(patches[patchMother].pos_of_first_ind + father);
            f = relationship[genCount%2][std::max(fathers.back(), mothers.back())][std::min(fathers.back(), mothers.back())];

        }

        juveniles[whr].emplace_back((patches[patchMother].population[mother_PosInPatch].s +
                                     patches[patchMother].population[father].s)/2,
                                    (patches[patchMother].population[mother_PosInPatch].d +
                                     patches[patchMother].population[father].d)/2, f);
    }
}

void World::mutation(Individual& IndToMutate)
{
    std::uniform_real_distribution<double> unif(0, 1);

    /* Y a-t-il mutation ? */
    if(unif(generator) < mu)
    {
        /* Pour «retenir» quel trait doit muter
        false: d     true: s */
        bool sWasChosen = false;
        double trait = IndToMutate.d;

        /* On choisit quel trait mute */
        if(unif(generator) >= d_s_relativeMutation)
        {
            trait = IndToMutate.s;
            sWasChosen = true;
        }

        switch(typeMut)
        {
            case gaussian:
                trait = gaussMutation(trait);
                break;

            case uniform:
                trait = unifMutation(trait);
                break;
        }

        if (sWasChosen) {IndToMutate.s = trait;}
        else {IndToMutate.d = trait;}
    }
}

double World::gaussMutation(double t)
{
    std::normal_distribution<double> gauss(0,sigmaZ);
    double deltaMu = gauss(generator);

    return t*exp(deltaMu)/(expm1(deltaMu)*t + 1); //expm1(x) renvoie exp(x) - 1.
}

double World::unifMutation(double t)
{
    double lowerBound = t - sigmaZ;
    double upperBound = t + sigmaZ;

    if(lowerBound < 0) {lowerBound = 0;}
    if(upperBound > 1) {upperBound = 1;}

    std::uniform_real_distribution<double> unif(lowerBound, upperBound);

    return unif(generator);
}

int World::getFather(int patchMother, int mother)
{
    int father = 0;

    std::uniform_int_distribution<int> unif(0, patches[patchMother].K-1);

    do {father = unif(generator);}
    while (father == mother); //Pas de pseudo allofécondation

    return father;
}

bool World::redefinePollination(double p)
{
    std::uniform_real_distribution<double> unif(0, 1);

    if (unif(generator) <= p)
    {
        return true;
    }

    return false;
}

void World::calcNewRelationships(void)
{
    int i = 0, j = 0;

    for(i=0; i<int(globalPop.size()); i++)
    {
        for(j=0; j<=i; j++)
        {
            if (i == j)
            {
                /* On a besoin du taux de consanguinité de l'individu. */
                relationship[(genCount+1)%2][i][j] = 1/2 + patches[globalPop[i].patch].population[globalPop[i].posInPatch].f/2;
            }

            relationship[(genCount+1)%2][i][j] = (relationship[genCount%2][mothers[i]][mothers[j]] +
                                                  relationship[genCount%2][fathers[i]][fathers[j]] +
                                                  relationship[genCount%2][fathers[i]][mothers[j]] +
                                                  relationship[genCount%2][mothers[i]][fathers[j]])/4;
        }
    }

    /*On vide les vecteurs. */
    fathers.clear();
    mothers.clear();
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
    report << "Gestion de l'apparentement:" << relationshipIsManaged << std::endl;
    report << "Delta=" << delta << " c=" << c << std::endl;
    report << "Loi pour la mutation:" << typeMut << " mu=" << mu << " sigmaZ=" << sigmaZ;
    report << " Taux de mutationt relatif d/s=" << d_s_relativeMutation << std::endl;
    report << "Kmin=" << Kmin << " Kmax=" << Kmax << " SigmaK=" << sigmaK << std::endl;
    report << "Pmin=" << Pmin << " Pmax=" << Pmax << " SigmaP=" << sigmaP << std::endl;
    report << "Gen\tPatch\tInd\ts\td" << std::endl;
}

void World::writeReport(void)
{
    int i = 0, j = 0;


    for(j=0; j<NPatch; j++)
    {
        for(i=0; i<patches[j].K; i++)
        {
            report << genCount << '\t';
            report << j << '\t';
            report << i << '\t';
            report << patches[j].population[i].s << '\t';
            report << patches[j].population[i].d << std::endl;
        }
    }


    if (NGen - genCount < genReport) {report.close();}
}

