#include <vector>
#include <iostream>
#include <fstream>

#include "engine.h"
#include "world.h"

int getParamsFromTxt(std::vector<double>& params)
{
    int NWorld, i;
    std::ifstream config("config.txt");

    config.ignore(500,'>');
    config.ignore(256,' ');
    config >> NWorld;

    for (i=0; i<16*NWorld; i++)
    {
        double val;
        config.ignore(256,' ');
        config >> val;
        params.push_back(val);
    }

    return NWorld;
}

void runSimu(void)
{
    int i = 0;
    std::vector<double> params;

    int NWorld = getParamsFromTxt(params);

    for (i=0; i<NWorld; i++)
    {
        World world(i, params[i], params[NWorld+i], params[2*NWorld+i], params[3*NWorld+i], params[4*NWorld+i],
                    params[5*NWorld+i], params[6*NWorld+i], params[7*NWorld+i], params[8*NWorld+i],
                    params[9*NWorld+i], params[10*NWorld+i], params[11*NWorld+i], params[12*NWorld+i],
                    params[13*NWorld+i], params[14*NWorld+i], params[15*NWorld+i]);
        world.run();
    }
}
