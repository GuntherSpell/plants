#include <random>
#include <vector>
#include <chrono>

#include "patch.h"
#include "individual.h"

Patch::Patch(double p, int K, double sInit, double dInit, int pos_of_first_ind)
{
    int i = 0;

    this->p = p;
    this->K = K;

    this->pos_of_first_ind = pos_of_first_ind;

    /* La ligne suivante permet de réserver
    de la mémoire pour éviter les réallocations
    qui peuvent diminuer les performances. */
    population.reserve(K);

    for(i=0; i<K; i++)
    {
        /* On initialise des individus non consanguins. */
        population.emplace_back(sInit, dInit, 0);
    }

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed (seed);

    isPollenized();
}

void Patch::isPollenized(void)
{
    std::uniform_real_distribution<double> unif(0, 1);

    pollenized = false;
    if (unif(generator) <= p)
    {
        pollenized = true;
    }
}

void Patch::getPression(double delta, double c, bool dispNeeded, std::vector<double>& press)
{
    int i = 0;

    if (dispNeeded)
    {
        if(dispSeeds.empty()) //Si le vecteur est vide, il faut le remplir.
        {
            dispSeeds.reserve(2*K);
            for(i=0; i<K; i++)
            {
                population[i].calcPress(delta, c, pollenized, dispNeeded, dispSeeds);
            }
        }

        press.insert(press.end(), dispSeeds.begin(), dispSeeds.end());
    }

    else
    {
        for(i=0; i<K; i++)
        {
            population[i].calcPress(delta, c, pollenized, dispNeeded, press);
        }
    }
}
