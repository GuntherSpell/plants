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

int Patch::check_convergence(int reportCount)
{
    int i = 0;

    std::vector<double> trait_values;
    std::array<double, 2> new_mean_var = {0,0};

    trait_values.reserve(K);

    if(reportCount<2)
    {
        /* Pour d */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].d);
        }

        calc_mean_var(trait_values, new_mean_var);
        previous_d_means[reportCount%2] = new_mean_var[0];

        trait_values.clear();

        /* Pour s */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].s);
        }

        calc_mean_var(trait_values, new_mean_var);
        previous_s_means[reportCount%2] = new_mean_var[0];
    }

    else
    {
        if(d_hasConverged && s_hasConverged)
        {
            return 1; //Le patch a convergé, on renvoie 1 pour sommer tous le patchs convergé.
        }

        /* Pour d */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].d);
        }

        calc_mean_var(trait_values, new_mean_var);

        d_hasConverged = check_stats(previous_d_means, previous_d_vars, new_mean_var);

        previous_d_means[reportCount%2] = new_mean_var[0];

        trait_values.clear();

        /* Pour s */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].s);
        }

        calc_mean_var(trait_values, new_mean_var);

        s_hasConverged = check_stats(previous_s_means, previous_s_vars, new_mean_var);

        previous_s_means[reportCount%2] = new_mean_var[0];
    }

    return 0;
}

bool Patch::check_stats(std::array<double, 2> previous_means, std::array<double, 2> previous_vars, std::array<double, 2> new_vals)
{
    if(std::abs(new_vals[0] - previous_means[0])/previous_means[0] < 0.05 &&
       std::abs(new_vals[0] - previous_means[1])/previous_means[1] < 0.05)
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
    /*
    double var = 0;

    for(i=0; i<int(values.size()); i++)
    {
        var += (values[i] - mean)*(values[i] - mean);
    }
    */
    mean_var[0] = mean;
    //mean_var[1] = var;
}
