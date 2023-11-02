#include "setup.hpp"
#include "FiniteDifferenceMethod.hpp"

double foo(double x, double y){
    return 1;
};
double fooN(double x, double y){
    return 0;
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
    FiniteDifferenceMethod drv(N);
    // drv.s = 1;
    // drv.create_blueprint();
    // drv.Source(f);
    // drv.fast_boundary(gD,gN);
    // drv.read_blueprint();
    drv.solve(f,gD,gN,true);

    // cout << drv.A << endl;
    // cout << drv.rhs <<endl;
    cout << drv.result.maxCoeff() << endl;
    return 0;
}
