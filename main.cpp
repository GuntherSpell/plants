#include <sstream>
#include <array>

#include "world.h"

int main(int argc, char *argv[])
{
    int i = 0;
    std::array<double, 29> params;

    int checkSum = 0; // égal à la longeur de params si les conversions ont bien eu lieu.

    /* Conversion char* -> double */
    for(i=0; i<29; i++)
    {
        std::istringstream iss(argv[i+1]);

        if(iss >> params[i])
        {
            checkSum ++;
        }
    }


    if(checkSum == params.size())
    {
        World world(params[0], params[1], params[2], params[3],
        params[4], params[5], params[6], params[7], params[8], params[9],
        params[10], params[11], params[12], params[13], params[14],
        params[15], params[16], params[17], params[18], params[19],
        params[20], params[21], params[22], params[23], params[24], params[25],
        params[26], params[27], params[28]);
        world.run(params[0]);
    }

    return 0;
}
