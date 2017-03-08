#include <random>
#include <cmath>
#include <chrono>

#include "individual.h"

Individual::Individual(double s, double d)
{
    this->s = s;
    this->d = d;
}

void Individual::mutation(double mu, double sigmaZ)
{
    double deltaMu = 0;

    /* Seed basée sur l'horloge pour le générateur de nombres aléatoires */
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();

    /* Générateur de nombres pseudo-aléatoires */
    std::default_random_engine generator (seed);

    /* On crée une série uniforme entre 0 et 1 */
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    /* Y a-t-il mutation ? */
    if(unif(generator) < mu)
    {
        /* Si oui, on choisit aléatoirement le degré de mutation
        Pour cela, on crée une distribution normale en fonction de sigma Z */
        std::normal_distribution<double> gauss(0,sigmaZ);
        deltaMu = gauss(generator);

        /* On choisit quel trait mute */
        if(unif(generator) > 0.5)
        {
            s = s*exp(deltaMu)/((exp(deltaMu) - 1)*s + 1);
        }
        else
        {
            d = d*exp(deltaMu)/((exp(deltaMu) - 1)*d + 1);
        }
    }
}

void Individual::calcPress(double delta, double c, bool pollenized, bool disp, std::vector<double>& press)
{
    if(disp)
    {
        press.push_back(s*(1-delta)*(d/2)*(1-c));

        if(pollenized) {press.push_back((1-s)*(1-c)*(d/2));}

        else {press.push_back(0);}
    }

    else
    {
        press.push_back(s*(1-delta)*(1-d));

        if(pollenized) {press.push_back((1-s)*(1-d));}

        else {press.push_back(0);}
    }
}

