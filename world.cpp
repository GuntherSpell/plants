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

World::World(int idWorld, int NPatch, double delta, double c, bool relatednessIsManaged, double mitigateRelatedness,
             bool rangeToBeShifted, int shiftFrequency,
             int typeMut, double mu, double sigmaZ, double d_s_relativeMutation, int Kdistr, int Kmin, int Kmax, int sigmaK,
             int Pdistr, double Pmin, double Pmax, double sigmaP, double sInit, double dInit,
             bool convergenceToBeChecked, int NPatchToConverge, int NGenToConverge, double relativeConvergence,
             double absoluteConvergence, int checkConvergenceFrequency, int NGen, int genReport, bool logPoll_is_to_be_written)
{
    int i = 0, j = 0;

    this->NPatch = NPatch;

    this->delta = delta;
    this->c = c;

    this->relatednessIsManaged = relatednessIsManaged;
    this->mitigateRelatedness = mitigateRelatedness;

    this->rangeToBeShifted = rangeToBeShifted;
    this->shiftFrequency = shiftFrequency;

    this->typeMut = distrMut(typeMut);
    this->mu = mu;
    this->sigmaZ = sigmaZ;
    this->d_s_relativeMutation = d_s_relativeMutation;

    this->convergenceToBeChecked = convergenceToBeChecked;
    this->NPatchToConverge = NPatchToConverge;
    this->NGenToConverge = NGenToConverge;
    this->relativeConvergence = relativeConvergence;
    this->absoluteConvergence = absoluteConvergence;
    this->checkConvergenceFrequency = checkConvergenceFrequency;

    this->NGen = NGen;
    genCount = 0;

    n_choose_2 = 0;

    report.open ("report_" + std::to_string(idWorld) + ".txt");
    this->genReport = genReport;

    this->logPoll_is_to_be_written = logPoll_is_to_be_written;

    if(logPoll_is_to_be_written)
    {
        logPoll.open("logPoll_" + std::to_string(idWorld) + ".txt");
    }

    if(relatednessIsManaged)
    {
        relation_report.open("relation_" + std::to_string(idWorld) + ".txt");
    }

    /* Les lignes suivantes permettent de réserver
    de la mémoire pour éviter les réallocations
    qui peuvent diminuer les performances. */
    patches.reserve(NPatch);
    juveniles[0].reserve(Kmax);
    juveniles[1].reserve(Kmax);

    /* Deux vecteurs qui permettent de stocker les valeurs de K et P
    pour chaque patch avant de construire les patchs. */
    std::vector<int> list_of_K;
    std::vector<double> list_of_P;
    list_of_K.reserve(NPatch);
    list_of_P.reserve(NPatch);

    if(Kdistr == 0)
    {
        for(i=0; i<NPatch; i++)
        {
            list_of_K.push_back(GaussDistr(Kmin, Kmax, sigmaK, i));
        }
    }
    else
    {
        int K_to_reach = 0; // La somme des K si on avait une distribution gaussienne.

        for(i=0; i<NPatch; i++)
        {
            K_to_reach += GaussDistr(Kmin, Kmax, sigmaK, i);
        }

        double Kstep = (K_to_reach - Kmin*NPatch)/((NPatch - 1)*(0.5 + (NPatch - 1)/2));
        /* ((NPatch - 1)*(0.5 + (NPatch - 1)/2)) est la somme des entiers de 1 à NPatch - 1
        Ce nombre correspond aux nombres de fois qu'on doit ajouter Kstep à la population. */

        if(Kdistr == 1)
        {
            for(i=0; i<NPatch; i++)
            {
                list_of_K.push_back(Kmin + i*Kstep);
            }
        }

        else
        {
            for(i=0; i<NPatch; i++)
            {
                list_of_K.push_back(Kmin + (30-i)*Kstep);
            }
        }
    }

    if(Pdistr == 0)
    {
        for(i=0; i<NPatch; i++)
        {
            list_of_P.push_back(GaussDistr(Pmin, Pmax, sigmaP, i));
        }
    }
    else
    {
        double P_to_reach = 0; // La somme des P si on avait une distribution gaussienne.

        for(i=0; i<NPatch; i++)
        {
            P_to_reach += GaussDistr(Pmin, Pmax, sigmaP, i);
        }

        double Pstep = (P_to_reach - Pmin*NPatch)/((NPatch - 1)*(0.5 + (NPatch - 1)/2));
        /* ((NPatch - 1)*(0.5 + (NPatch - 1)/2)) est la somme des entiers de 1 à NPatch - 1
        Ce nombre correspond aux nombres de fois qu'on doit ajouter Pstep à la population. */

        if(Pdistr == 1)
        {
            for(i=0; i<NPatch; i++)
            {
                list_of_P.push_back(Pmin + i*Pstep);
            }
        }

        else
        {
            for(i=0; i<NPatch; i++)
            {
                list_of_P.push_back(Pmin + (30-i)*Pstep);
            }
        }
    }

    /* Construction des patchs */
    int Ktot = 0;
    double Ptot = 0;

    for(i=0; i<NPatch; i++)
    {
        patches.emplace_back(list_of_P[i], list_of_K[i], sInit, dInit, Ktot);
        Ktot += patches[i].K;
        Ptot += patches[i].p;
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

    /* Préparation des vecteurs des moyennes des traits pour vérifier la convergence */
    if(convergenceToBeChecked)
    {
        n_choose_2 = factorial(NGenToConverge)/(2*factorial(NGenToConverge - 2));
        for(i=0; i<NPatch; i++)
        {
            patches[i].previous_d_means.assign(NGenToConverge - 1, 0);
            patches[i].previous_s_means.assign(NGenToConverge - 1, 0);

            for(j=(NGenToConverge - 1); j>0; j--)
            {
                patches[i].prev_gens_s_similarity_matrix.emplace_back(j);
                patches[i].prev_gens_d_similarity_matrix.emplace_back(j);
            }
        }
    }

    /* Préparation des matrices d'apparentement si nécessaire */
    if(relatednessIsManaged)
    {
        relatedness[0].reserve(Ktot);
        relatedness[1].reserve(Ktot);
        fathers.reserve(Ktot);
        mothers.reserve(Ktot);

        /* On part d'individus non apparentés. */
        for(i=0; i<Ktot; i++)
        {
            /* Pour chacune des deux matrices, on construit Ktot lignes remplies de i+1 zéros. */
            relatedness[0].emplace_back(i + 1);
            relatedness[1].emplace_back(i + 1);

            /* Pour la première matrice, il faut remplir la diagonale de 0.5. */
            relatedness[0][i][i] = 0.5;
        }
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed (seed);

    writeHeaders(Kdistr, Kmin, Kmax, sigmaK, Ktot, Pdistr, Pmin, Pmax, sigmaP, Ptot);
}

double World::GaussDistr(double minVal, double maxVal, double sigma, int posPatch)
{
    return (minVal + ((maxVal - minVal) * exp( - ((posPatch-(NPatch/2))*(posPatch-(NPatch/2))) / (2*sigma*sigma))));
}

void World::run(int idWorld)
{
    int i = 0, progress = 0, checkCount = 0;

    std::cout << "Progression du monde " << idWorld << " :" << std::endl;
    printProgress(0);

    for(genCount=0; genCount<=NGen; genCount++)
    {
        if(genCount%genReport == 0)
        {
            writeReport();
        }

        for(i=0; i<NPatch; i++)
        {
            patches[i].pollenized = redefinePollination(patches[i].p);
        }

        if(logPoll_is_to_be_written)
        {
            writeLogPoll();
        }

        if(rangeToBeShifted && genCount%shiftFrequency == 0)
        {
            /* Déplacement des populations */
            for(i=0; i < NPatch - 1; i++)
            {
                patches[i].population = patches[i+1].population;

                /* On tue les individus en trop. */
                while(patches[i].K < int(patches[i].population.size()))
                {
                    patches[i].population.pop_back();
                }
            }

            /* Le patch de droite (le nouveau) est vidé. */
            patches[NPatch - 1].population.clear();

            /* Réévaluer la position du premier individu dans le patch (par rapport à la pop globale). */
            unsigned int new_pos_first_ind = patches[0].population.size();
            for(i=1; i<NPatch; i++)
            {
                patches[i].pos_of_first_ind = new_pos_first_ind;
                new_pos_first_ind += patches[i].population.size();
            }

            /* Il faut recréer la globalPop avant de pouvoir créer la nouvelle génération. */
            globalPop.clear();

            for(i=0; i<NPatch; i++)
            {
                for(int j = 0; j < int(patches[i].population.size()); j++)
                {
                    IndividualPosition InfoToAdd;   //
                    InfoToAdd.patch = i;            // On ne peut pas ajouter l'info dans le vecteur sans la construire avant.
                    InfoToAdd.posInPatch = j;       //

                    globalPop.push_back(InfoToAdd);
                }
            }

            for(i=0; i<NPatch; i++)
            {
                createNextGen(i);
            }

            /* D'une fois que la nouvelle génération est créée, le monde a retrouvé sa population normale.
            Il faut donc recréer globalPop une nouvelle fois. */
            globalPop.clear();

            for(i=0; i<NPatch; i++)
            {
                for(int j = 0; j < int(patches[i].population.size()); j++)
                {
                    IndividualPosition InfoToAdd;   //
                    InfoToAdd.patch = i;            // On ne peut pas ajouter l'info dans le vecteur sans la construire avant.
                    InfoToAdd.posInPatch = j;       //

                    globalPop.push_back(InfoToAdd);
                }
            }
        }

        else
        {
            for(i=0; i<NPatch; i++)
            {
                createNextGen(i);
            }
        }


        if(convergenceToBeChecked && genCount%checkConvergenceFrequency == 0 &&
           (genCount >= 100000 - (NGenToConverge - 1)*checkConvergenceFrequency))
        {
            int sumOfConvergedPatches = 0;

            for(i=0; i<NPatch; i++)
            {
                sumOfConvergedPatches +=
                patches[i].check_convergence(checkCount, NGenToConverge, n_choose_2, relativeConvergence, absoluteConvergence);
            }

            if(sumOfConvergedPatches >= NPatchToConverge)
            {
                if(relatednessIsManaged)
                {
                    writeRelatednesses();
                }

                return; // Si on a rempli le critère de convergence, on arrête la simu.
            }

            checkCount ++;
        }

        if(relatednessIsManaged)
        {
            calcNewRelatednesses();
        }

        /* Indique la progression à l'écran */
        if (genCount*78/NGen > progress)
        {
            progress = genCount*78/NGen;
            printProgress(progress);
        }
    }

    if(relatednessIsManaged)
    {
        writeRelatednesses();
    }
}

void World::createNextGen(int idPatch)
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

    for(i=0; i<patches[idPatch].K; i++)
    {
        int chosenMother = weighted(generator);

        /* Une mère fait de l'autof si elle a la même parité que la première mère. */
        bool autof = false;
        if (chosenMother%2 == firstMother%2)
        {
            autof = true;
        }

        /* Pour reprendre l'id de la mère (on en avait le double car autof et allof). */
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

    /* On récupère la position relative de la mère dans son patch. */
    int patchMother = globalPop[mother].patch;
    int mother_PosInPatch = globalPop[mother].posInPatch;

    /* Issue d'autof */
    if(autof)
    {
        double f = 0.5;

        if(relatednessIsManaged)
        {
            f = 0.5 + patches[patchMother].population[mother_PosInPatch].f*0.5;
            mothers.push_back(mother);
            fathers.push_back(mother);
        }

        juveniles[whr].emplace_back(patches[patchMother].population[mother_PosInPatch].s,
                                    patches[patchMother].population[mother_PosInPatch].d, f);

    }

    /* Sinon, on cherche un père. */
    else
    {
        double f = 0;

        int father = getFather(patchMother, mother_PosInPatch);

        if(relatednessIsManaged)
        {
            mothers.push_back(mother);
            fathers.push_back(patches[patchMother].pos_of_first_ind + father);
            f = relatedness[genCount%2][std::max(fathers.back(), mothers.back())][std::min(fathers.back(), mothers.back())];

        }

        juveniles[whr].emplace_back(0.5*(patches[patchMother].population[mother_PosInPatch].s +
                                     patches[patchMother].population[father].s),
                                    0.5*(patches[patchMother].population[mother_PosInPatch].d +
                                     patches[patchMother].population[father].d), f);
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

        if (sWasChosen)
        {
            IndToMutate.s = trait;
        }
        else
        {
            IndToMutate.d = trait;
        }
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

    if(lowerBound < 0)
    {
        lowerBound = 0;
    }

    if(upperBound > 1)
    {
        upperBound = 1;
    }

    std::uniform_real_distribution<double> unif(lowerBound, upperBound);

    return unif(generator);
}

int World::getFather(int patchMother, int mother)
{
    int father = 0;

    std::uniform_int_distribution<int> unif(0, patches[patchMother].population.size() - 1);

    do
    {
        father = unif(generator);
    }
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

void World::calcNewRelatednesses(void)
{
    int i = 0, j = 0;

    /* Les apparentements sont multipliés par 1 - mu pour introduire le fait
       qu'une partie du génome change à cause des mutations. */
    for(i=0; i<int(globalPop.size()); i++)
    {
        for(j=0; j<=i; j++)
        {
            /* Il faut remplir la diagonale pour les indivdus ayant un ou deux parents en commun. */
            if(i == j)
            {
                /* On a besoin du taux de consanguinité de l'individu. */
                relatedness[(genCount+1)%2][i][j] = (1 - mitigateRelatedness) *
                (0.5 + 0.5*patches[globalPop[i].patch].population[globalPop[i].posInPatch].f);
            }

            else
            {
                relatedness[(genCount+1)%2][i][j] = (1 - mitigateRelatedness) *
                (relatedness[genCount%2][std::max(mothers[i], mothers[j])][std::min(mothers[i], mothers[j])] +
                 relatedness[genCount%2][std::max(fathers[i], fathers[j])][std::min(fathers[i], fathers[j])] +
                 relatedness[genCount%2][std::max(fathers[i], mothers[j])][std::min(fathers[i], mothers[j])] +
                 relatedness[genCount%2][std::max(mothers[i], fathers[j])][std::min(mothers[i], fathers[j])])*0.25;
            }
        }
    }

    /* On vide les vecteurs. */
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

void World::writeHeaders(int Kdistr, int Kmin, int Kmax, int sigmaK, int Ktot,
                         int Pdistr, double Pmin, double Pmax, double sigmaP, double Ptot)
{
    report << "Nombre de patchs=" << NPatch << std::endl;
    report << "Gestion de l'apparentement:" << relatednessIsManaged;
    report << " Correction apparentement=" << mitigateRelatedness << std::endl;
    report << "Delta=" << delta << " c=" << c << std::endl;
    report << "Loi pour la mutation:" << typeMut << " mu=" << mu << " sigmaZ=" << sigmaZ;
    report << " Taux de mutationt relatif d/s=" << d_s_relativeMutation << std::endl;
    report << "KDistr:" << Kdistr << " Kmin=" << Kmin << " Kmax=" << Kmax << " SigmaK=" << sigmaK << " K_tot=" << Ktot << std::endl;
    report << "PDistr:" << Pdistr << " Pmin=" << Pmin << " Pmax=" << Pmax << " SigmaP=" << sigmaP << " P_tot=" << Ptot << std::endl;
    report << "Vérifier convergence:" << convergenceToBeChecked << " N patchs à converger=" << NPatchToConverge;
    report << " N gen à converger=" << NGenToConverge << " Relatif=" << relativeConvergence << " Absolu=" << absoluteConvergence;
    report << " Fréquence=" << checkConvergenceFrequency << std::endl;
    report << "Shift=" << rangeToBeShifted << " Fréquence=" << shiftFrequency << std::endl;
    report << "Gen\tPatch\tInd\ts\td" << std::endl;

    if(logPoll_is_to_be_written)
    {
        logPoll << "Gen" << '\t' << "Patch" << '\t' << "Etat" << std::endl;
    }
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
            report << std::round(patches[j].population[i].s * 1000) / 1000 << '\t';
            report << std::round(patches[j].population[i].d * 1000) / 1000 << std::endl;
        }
    }


    if (NGen - genCount < genReport)
    {
        report.close();
    }
}

void World::writeLogPoll(void)
{
    int i = 0;

    for(i=0; i<NPatch; i++)
    {
        logPoll << genCount << '\t' << i << '\t' << patches[i].pollenized << std::endl;
    }
}

void World::writeRelatednesses(void)
{
    int i = 0, j = 0, k = 0;
    int ind_abs_id = 0; // La position absolue de l'individu (la ligne) concerné.

    relation_report << '\t';

    /* Première ligne. */
    for(i=0; i<NPatch; i++)
    {
        for(j=0; j<patches[i].K; j++)
        {
            relation_report << j << '\t';
        }
    }

    /* Reste de la matrice. */
    for(i=0; i<NPatch; i++)
    {
        for(j=0; j<patches[i].K; j++)
        {
            relation_report << std::endl << j << '\t';

            for(k=0; k<=ind_abs_id; k++)
            {
                relation_report << relatedness[0][ind_abs_id][k] << '\t';
            }

            ind_abs_id ++;
        }
    }
}

int World::factorial(int n)
{
    if(n <= 1)
    {
        return 1;
    }

    return n*factorial(n - 1);
}
