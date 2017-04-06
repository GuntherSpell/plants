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

    d_hasConverged = false;
    s_hasConverged = false;

    previous_d_means = {0,0};
    previous_s_means = {0,0};
}

void Patch::getDispPress(double delta, double c, std::vector<double>& press)
{
    if(dispSeeds.empty()) //Si le vecteur est vide, il faut le remplir.
    {
        int i = 0;

        dispSeeds.reserve(2*K);
        for(i=0; i<K; i++)
        {
            population[i].calcDispPress(delta, c, pollenized, dispSeeds);
        }
    }

    press.insert(press.end(), dispSeeds.begin(), dispSeeds.end());
}

void Patch::getResidPress(double delta, std::vector<double>& press)
{
    int i = 0;

    for(i=0; i<K; i++)
    {
        population[i].calcResidPress(delta, pollenized, press);
    }
}

int Patch::check_convergence(int reportCount, double relativeConvergence, double absoluteConvergence)
{
    int i = 0;

    std::vector<double> trait_values;

    trait_values.reserve(K);

    if(reportCount<2)
    {
        /* Pour d */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].d);
        }

        previous_d_means[reportCount] = calc_mean(trait_values);

        trait_values.clear();

        /* Pour s */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].s);
        }

        previous_s_means[reportCount] = calc_mean(trait_values);
    }

    else
    {
        /* Pour d */
        if(!d_hasConverged)
        {
            for(i=0; i<K; i++)
            {
                trait_values.push_back(population[i].d);
            }

            double new_d_mean = calc_mean(trait_values);

            d_hasConverged = check_stats(previous_d_means, new_d_mean, relativeConvergence, absoluteConvergence);

            previous_d_means[reportCount%2] = new_d_mean;

            trait_values.clear();
        }

        /* Pour s */
        if(!s_hasConverged)
        {
            for(i=0; i<K; i++)
            {
                trait_values.push_back(population[i].s);
            }

            double new_s_mean = calc_mean(trait_values);

            s_hasConverged = check_stats(previous_s_means, new_s_mean, relativeConvergence, absoluteConvergence);

            previous_s_means[reportCount%2] = new_s_mean;
        }

        if(d_hasConverged && s_hasConverged)
        {
            return 1; //Le patch a convergé, on renvoie 1 pour sommer tous le patchs convergés.
        }
    }

    return 0;
}

bool Patch::check_stats(std::array<double, 2> previous_means, double new_mean, double relativeConvergence, double absoluteConvergence)
{
    /* On a deux critères avec un OU logique
    En variation relative: efficace quand la valeur est haute
    En variation absolue: efficace quand la valeur est basse */
    if((std::abs(new_mean - previous_means[0])/previous_means[0] < relativeConvergence &&
       std::abs(new_mean - previous_means[1])/previous_means[1] < relativeConvergence) ||
       (std::abs(new_mean - previous_means[0]) < absoluteConvergence &&
       std::abs(new_mean - previous_means[1]) < absoluteConvergence))
    {
        return true;
    }

    return false;
}

double Patch::calc_mean(std::vector<double> values)
{
    int i = 0;

    double sum = 0;

    for(i=0; i<int(values.size()); i++)
    {
        sum += values[i];
    }

    return sum/values.size();
}
