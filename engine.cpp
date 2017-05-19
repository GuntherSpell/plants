#include <vector>
#include <array>
#include <fstream>

#include "engine.h"
#include "world.h"

void getParamsFromTxt(std::vector<std::array<double, 25>>& params, int& NWorld, int& NReplicats)
{
    int i = 0, j = 0;
    std::ifstream config("config.txt");

    config.ignore(256,'\n');
    config >> NWorld;

    config.ignore(256,'\t');
    config >> NReplicats;

    params.reserve(NWorld);

    config.ignore(256,'\n');
    config.ignore(256,'\n');

    for(i=0; i<NWorld; i++)
    {
        std::array<double, 25> paramOneWorld = {0,};

        for(j=0; j<25; j++)
        {
            double val = 0;
            config >> val;
            paramOneWorld[j] = val;
            config.ignore(1,'\t');
        }

        params.push_back(paramOneWorld);
    }
}

void runSimu(void)
{
    int i = 0, j = 0, NWorld = 0, NReplicats = 0;
    std::vector<std::array<double, 25>> params;

    getParamsFromTxt(params, NWorld, NReplicats);

    for(i=0; i<NWorld; i++)
    {
        for(j=0; j<NReplicats; j++)
        {
            World world(i*NReplicats + j, params[i][0], params[i][1], params[i][2], params[i][3],
            params[i][4], params[i][5], params[i][6], params[i][7], params[i][8], params[i][9],
            params[i][10], params[i][11], params[i][12], params[i][13], params[i][14],
            params[i][15], params[i][16], params[i][17], params[i][18], params[i][19],
            params[i][20], params[i][21], params[i][22], params[i][23], params[i][24]);
            world.run(i*NReplicats + j);
        }
    }
}
