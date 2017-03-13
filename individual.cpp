#include <random>
#include <cmath>
#include <chrono>

#include "individual.h"

Individual::Individual(double s, double d)
{
    this->s = s;
    this->d = d;

    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    generator.seed (seed);
}

void Individual::mutation(double mu, double sigmaZ, distrMut typeMut)
{
    /* Pour «retenir» quel trait doit muter
       false: d     true: s */
    bool sWasChosen = false;
    double t = d;

    /* On crée une série uniforme entre 0 et 1 */
    std::uniform_real_distribution<double> unif(0.0, 1.0);

    /* Y a-t-il mutation ? */
    if(unif(generator) < mu)
    {
        /* On choisit quel trait mute */
        if(unif(generator) >= 0.5)
        {
            t = s;
            sWasChosen = true;
        }

        switch(typeMut)
        {
            case gaussian:
                t = gaussMutation(sigmaZ, t);
                break;

            case uniform:
                t = unifMutation(sigmaZ, t);
                break;
        }

        if (sWasChosen) {s = t;}
        else {d = t;}
    }
}

double Individual::gaussMutation (double sigmaZ, double t)
{
    std::normal_distribution<double> gauss(0,sigmaZ);
    double deltaMu = gauss(generator);

    return t*exp(deltaMu)/((exp(deltaMu) - 1)*t + 1);
}

double Individual::unifMutation (double sigmaZ, double t)
{
    double lowerBound = t - sigmaZ;
    double upperBound = t + sigmaZ;

    if(lowerBound < 0) {lowerBound = 0;}
    if(upperBound > 1) {upperBound = 1;}

    std::uniform_real_distribution<double> unif(lowerBound, upperBound);

    return unif(generator);
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

