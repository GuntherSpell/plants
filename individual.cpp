#include <random>
#include <cmath>
#include <chrono>

#include "individual.h"

Individual::Individual(double s, double d)
{
    this->s = s;
    this->d = d;
}

void Individual::calcPress(double delta, double c, bool pollenized, bool dispNeeded, std::vector<double>& press)
{
    if(dispNeeded)
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

