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

        dispSeeds.reserve(2*population.size());
        for(i=0; i < population.size(); i++)
        {
            population[i].calcDispPress(delta, c, pollenized, dispSeeds);
        }
    }

    press.insert(press.end(), dispSeeds.begin(), dispSeeds.end());
}

void Patch::getResidPress(double delta, std::vector<double>& press)
{
    int i = 0;

    for(i=0; i < population.size(); i++)
    {
        population[i].calcResidPress(delta, pollenized, press);
    }
}

int Patch::check_convergence(int checkCount, int NGenToConverge, int n_choose_2, double relativeConvergence, double absoluteConvergence)
{
    int i = 0;

    std::vector<double> trait_values;

    trait_values.reserve(K);

    /* Avant d'avoir suffisamment de générations de référence, on remplit la matrice. */
    if(checkCount < NGenToConverge)
    {
        /* Pour d */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].d);
        }

        double new_d_mean = calc_mean(trait_values);

        for(i=0; i<checkCount; i++)
        {
            prev_gens_d_similarity_matrix[i][NGenToConverge - checkCount - 1] =
            check_stats(previous_d_means[i], new_d_mean, relativeConvergence, absoluteConvergence);
        }

        previous_d_means[checkCount] = new_d_mean;

        trait_values.clear();

        /* Pour s */
        for(i=0; i<K; i++)
        {
            trait_values.push_back(population[i].s);
        }

        double new_s_mean = calc_mean(trait_values);

        for(i=0; i<checkCount; i++)
        {
            prev_gens_s_similarity_matrix[i][NGenToConverge - checkCount - 1] =
            check_stats(previous_s_means[i], new_s_mean, relativeConvergence, absoluteConvergence);
        }

        previous_s_means[checkCount] = new_s_mean;
    }

    /* D'une fois que la matrice est remplie, on peut vérifier la convergence. */
    else
    {
        int j = 0;
        int convergeSum = 0;

        /* Pour d */
        if(!d_hasConverged)
        {
            /* On vérifie que ce trait a convergé. */
            for(i=0; i<(NGenToConverge - 1); i++)
            {
                for(j=0; j<(NGenToConverge - 1 - i); j++)
                {
                    convergeSum += prev_gens_d_similarity_matrix[i][j];
                }
            }

            if(convergeSum == n_choose_2)
            {
                d_hasConverged = 1;
            }

            else
            {
                for(i=0; i<K; i++)
                {
                    trait_values.push_back(population[i].d);
                }

                double new_d_mean = calc_mean(trait_values);

                while(i < checkCount%NGenToConverge)
                {
                    prev_gens_d_similarity_matrix[i][NGenToConverge - 1 - checkCount%NGenToConverge] =
                    check_stats(previous_d_means[i], new_d_mean, relativeConvergence, absoluteConvergence);
                    i++;
                }

                for(i=0; i < NGenToConverge - 1 - checkCount%NGenToConverge; i++)
                {
                    prev_gens_d_similarity_matrix[checkCount%NGenToConverge][i] =
                    check_stats(previous_d_means[NGenToConverge - 1 - i], new_d_mean, relativeConvergence, absoluteConvergence);
                }

                previous_d_means[checkCount%NGenToConverge] = new_d_mean;

                trait_values.clear();

            }
        }

        /* Pour s */
        if(!s_hasConverged)
        {
            /* On vérifie que ce trait a convergé. */
            convergeSum = 0;

            for(i=0; i<(NGenToConverge - 1); i++)
            {
                for(j=0; j<(NGenToConverge - 1 - i); j++)
                {
                    convergeSum += prev_gens_s_similarity_matrix[i][j];
                }
            }

            if(convergeSum == n_choose_2)
            {
                s_hasConverged = 1;
            }

            else
            {
                for(i=0; i<K; i++)
                {
                    trait_values.push_back(population[i].s);
                }

                double new_s_mean = calc_mean(trait_values);

                while(i < checkCount%NGenToConverge)
                {
                    prev_gens_s_similarity_matrix[i][NGenToConverge - 1 - checkCount%NGenToConverge] =
                    check_stats(previous_s_means[i], new_s_mean, relativeConvergence, absoluteConvergence);
                    i++;
                }

                for(i=0; i < NGenToConverge - 1 - checkCount%NGenToConverge; i++)
                {
                    prev_gens_s_similarity_matrix[checkCount%NGenToConverge][i] =
                    check_stats(previous_s_means[NGenToConverge - 1 - i], new_s_mean, relativeConvergence, absoluteConvergence);
                }

                previous_s_means[checkCount%NGenToConverge] = new_s_mean;
            }
        }

        if(d_hasConverged && s_hasConverged)
        {
            return 1; //Le patch a convergé, on renvoie 1 pour sommer tous les patchs convergés.
        }
    }

    return 0;
}

bool Patch::check_stats(double first_mean, double second_mean, double relativeConvergence, double absoluteConvergence)
{
    /* On a deux critères avec un OU logique
    En variation relative: efficace quand la valeur est haute
    En variation absolue: efficace quand la valeur est basse */
    if(std::abs(first_mean - second_mean)/second_mean < relativeConvergence ||
       (std::abs(first_mean - second_mean) < absoluteConvergence))
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
