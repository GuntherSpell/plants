#include <random>
#include <cmath>
#include <chrono>

#include "individual.h"

Individual::Individual(double s, double d, double f)
{
    this->s = s;
    this->d = d;

    this->f = f;
}

void Individual::calcDispPress(double delta, double c, bool pollenized, std::vector<double>& press)
{
    press.push_back(s*(1-delta)*(d/2)*(1-c));

    if(pollenized)
    {
        press.push_back((1-s)*(1-c)*(d/2));
    }

    else
    {
        press.push_back(0);
    }
}

void Individual::calcResidPress(double delta, bool pollenized, std::vector<double>& press)
{
    press.push_back(s*(1-delta)*(1-d));

    if(pollenized)
    {
        press.push_back((1-s)*(1-d));
    }

    else
    {
        press.push_back(0);
    }
}
