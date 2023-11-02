#include "FiniteElementMethod.hpp"

double func(double x, double y){
    return x + y*y;
}

int main(int argc, char const *argv[])
{
    int N = 5;
    double f(0),gD(0),gN(0);
    if(argc>1){
        N = atoi(argv[1]);
        if (argc>2){
            f = atof(argv[2]);
            if (argc>4){
                gD = atof(argv[3]);
                gN = atof(argv[4]);
            }
        }
    }
    
    ExoticFiniteElementMethod driver(N);
    driver.solve(f,gD,gN);
   

    return 0;
}
