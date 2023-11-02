#include "setup.hpp"
#include "FiniteDifferenceMethod.hpp"
#include "FiniteElementMethod.hpp"
#include "Convergence.hpp"


int main(int argc, char const *argv[])
{
    int N(5);
    double f(1);
    if(argc>1){
        N = atoi(argv[1]);
        if(argc>2){
            f = atof(argv[2]);
        }
    }
    Convergence analysis(5);
    analysis.compute_fdm(f,0,0);
    analysis.compute_norms();
    return 0;
}
