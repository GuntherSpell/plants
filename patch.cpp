#include <random>
#include <vector>
#include <chrono>

#include "patch.h"
#include "individual.h"

Patch::Patch(double p, int K, double sInit, double dInit)
{
    int i = 0;

    this->p = p;
    this->K = K;
    isPollenized();

    for(i=0; i<K; i++)
    {
        population.emplace_back(sInit, dInit);
    }
}

void Patch::isPollenized(void)
{
    /* Seed basée sur l'horloge pour le générateur de nombres aléatoires */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    /* Générateur de nombres pseudo-aléatoires */
    std::default_random_engine generator (seed);

    /* On crée une série uniforme entre 0 et 1 */
    std::uniform_real_distribution<double> unif(0, 1);

    if (unif(generator) <= p) {pollenized = true;}
    else {pollenized = false;}
}

void Patch::getPression(double delta, double c, bool disp, std::vector<double>& press)
{
    int i = 0;

    for(i=0; i<K; i++)
    {
        population[i].calcPress(delta, c, pollenized, disp, press);
    }
}
