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
    previous_d_vars = {0,0};
    previous_s_means = {0,0};
    previous_s_vars = {0,0};
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

bool Patch::check_stats(std::array<double, 2> previous_means, std::array<double, 2> previous_vars, std::array<double, 2> new_vals)
{
    if(std::abs(new_vals[0] - previous_means[0])/previous_means[0] < 0.015 &&
       std::abs(new_vals[0] - previous_means[1])/previous_means[1] < 0.015)
    {
        return true;
    }

    return false;
}

void Patch::calc_mean_var(std::vector<double> values, std::array<double, 2>& mean_var)
{
    int i = 0;

    /* Calcul de la moyenne */
    double sum = 0;

    for(i=0; i<int(values.size()); i++)
    {
        sum += values[i];
    }

    double mean = sum/values.size();

    /* Calcul de la variance */
    double var = 0;

    for(i=0; i<int(values.size()); i++)
    {
        var += (values[i] - mean)*(values[i] - mean);
    }

    mean_var[0] = mean;
    mean_var[1] = var;
}
