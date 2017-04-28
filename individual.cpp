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
    /* La dépression de consanguinité subie par l'individu, déterminée grace à la fonction f_to_delta. */
    double ind_delta = f_to_delta(delta, f);

    press.push_back(s*(1-ind_delta)*(1-c)*(d/2));

    if(pollenized)
    {
        press.push_back((1-s)*(1-ind_delta)*(1-c)*(d/2));
    }

    else
    {
        press.push_back(0);
    }
}

void Individual::calcResidPress(double delta, bool pollenized, std::vector<double>& press)
{
    /* La dépression de consanguinité subie par l'individu, déterminée grace à la fonction f_to_delta. */
    double ind_delta = f_to_delta(delta, f);

    press.push_back(s*(1-ind_delta)*(1-d));

    if(pollenized)
    {
        press.push_back((1-s)*(1-ind_delta)*(1-d));
    }

    else
    {
        press.push_back(0);
    }
}

double Individual::f_to_delta(double delta, double f)
{
    /* La différence entre la dépression consanguine pour f = 0.5 et f = 1. */
    double epsilon = 0.1;

    /* Attention à ne pas mettre delta <= à epsilon. */
    return (delta*f)/((epsilon/(delta - epsilon)) - (epsilon/(delta - epsilon) - 1)*f);
}
