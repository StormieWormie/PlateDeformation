#ifndef FDM_HPP
#define FDM_HPP

#include "setup.hpp"
#include "cluster.hpp"

/*
Solve:
c∇^4 u = f in Ω
u = g1 on Γ
δu/δn = g2 on Γ

stencil used:
      NN               1            r
   NW  N NE         2 -8  2       p n p
WW  W  C  E EE = 1 -8 20 -8 1 = r n m n r
   SW  S SE         2 -8  2       p n p
      SS               1            r
N=north, E=east, S=south, W=west, C=Center
constants m,n,p and r are for readability (q is skipped because my little brain cannot distinguish it from p)
*/

class FiniteDifferenceMethod
{
public:
//Attributes
    simple_cluster node_grid;
    int length, internal_length;
    int N_total, N_internal, N_boundary;
    double stepsize, s;
    vector<triplet> matrix_blueprint;
    SparseMatrix<double> A; VectorXd rhs;
    VectorXd result;

//Constructors
    FiniteDifferenceMethod(int N_=5):
    node_grid(N_),
    length(N_), internal_length(length-2),
    N_total(N_*N_), N_internal((N_-2)*(N_-2)), N_boundary(N_total-N_boundary),
    stepsize(1./double(N_-1)), s(1./stepsize),
    A(N_internal,N_internal),rhs(N_internal),
    result(N_total)
    {
        rhs.setZero();
        result.setZero();
    }

//Routines
    void create_blueprint(){
        double h = s*s*s*s;
        int I,J;
        double m(20*h), n(-8*h), p(2*h), r(1*h);
        //double m(4*h), n(3*h), p(2*h), r(1*h);
        for (int i = 0; i < N_internal; i++){
            matrix_blueprint.push_back(triplet(i,i,m));
        }
        for (int i = 0; i < internal_length; i++){
            for (int j = 0; j < internal_length-1; j++){
                I = i*internal_length+j;
                J = I+1;
                matrix_blueprint.push_back(triplet(I,J,n));
                matrix_blueprint.push_back(triplet(J,I,n));
            }
        }
        for (int i = 0; i < internal_length-1; i++){
            for (int j = 0; j < internal_length; j++){
                I = i*internal_length+j;
                J = I+internal_length;
                matrix_blueprint.push_back(triplet(I,J,n));
                matrix_blueprint.push_back(triplet(J,I,n));
            }
        }
        for (int i = 0; i < internal_length-1; i++){
            for (int j = 0; j < internal_length-1; j++){
                I = i*internal_length+j;
                J = I+internal_length+1;
                matrix_blueprint.push_back(triplet(I,J,p));
                matrix_blueprint.push_back(triplet(J,I,p));
            }
        }
        for (int i = 0; i < internal_length-1; i++){
            for (int j = 1; j < internal_length; j++){
                I = i*internal_length+j;
                J = I+internal_length-1;
                matrix_blueprint.push_back(triplet(I,J,p));
                matrix_blueprint.push_back(triplet(J,I,p));
            }
        }
        for (int i = 0; i < internal_length; i++){
            for (int j = 0; j < internal_length-2; j++){
                I = i*internal_length+j;
                J = I + 2;
                matrix_blueprint.push_back(triplet(I,J,r));
                matrix_blueprint.push_back(triplet(J,I,r));
            }
        }
        for (int i = 0; i < internal_length-2; i++){
            for (int j = 0; j < internal_length; j++){
                I = i*internal_length+j;
                J = I + 2*internal_length;
                matrix_blueprint.push_back(triplet(I,J,r));
                matrix_blueprint.push_back(triplet(J,I,r));
            }
        }
    }

    void Source(double value = 0){
        for (int i = 0; i < N_internal; i++){
            rhs(i) += value;
        }
    }

    void fast_boundary(double Dirichlet, double Neumann){
        double h = s*s*s*s;
        double m(20*h), n(-8*h), p(2*h), r(1*h);
        vector<int> ilist;
        for (int i = 0; i < internal_length; i++){
            ilist = {   internal_length+i,
                        (i+1)*internal_length-2,
                        N_internal-internal_length-1-i,
                        i*internal_length+1};
            for (int j : ilist){
                rhs(j) -= r*Dirichlet;
            }
        }
        
        for (int i = 0; i < internal_length; i++){
            ilist = {   i,
                        (i+1)*internal_length-1,
                        N_internal-1-i,
                        i*internal_length};
            for (int j : ilist){
                rhs(j) -= n*Dirichlet;
                rhs(j) -= r*(Neumann*stepsize + Dirichlet);
                rhs(j) -= 2*p*Dirichlet;
            }
        }
        ilist = {   0,
                    internal_length-1,
                    N_internal-internal_length,
                    N_internal-1};
        for (int i : ilist){
            rhs(i) += p*Dirichlet;
        }
    }

    void accurate_boundary(double Dirichlet, double Neumann){
        double h = s*s*s*s;
        double m(20*h), n(-8*h), p(2*h), r(1*h);
        vector<int> ilist;
        for (int i = 0; i < internal_length; i++){
            ilist = {   internal_length+i,
                        (i+1)*internal_length-2,
                        N_internal-internal_length-1-i,
                        i*internal_length+1};
            for (int j : ilist){
                rhs(j) -= r*Dirichlet;
            }
        }
        
        for (int i = 0; i < internal_length; i++){
            ilist = {   i,
                        (i+1)*internal_length-1,
                        N_internal-1-i,
                        i*internal_length};
            for (int j : ilist){
                rhs(j) -= n*Dirichlet;
                rhs(j) -= r*2*stepsize*Neumann;
                matrix_blueprint.push_back(triplet(j,j,r));
                rhs(j) -= 2*p*Dirichlet;
            }
        }
        ilist = {   0,
                    internal_length-1,
                    N_internal-internal_length,
                    N_internal-1};
        for (int i : ilist){
            rhs(i) += p*Dirichlet;
        }
    }
    
    void read_blueprint(){
        A.setFromTriplets(matrix_blueprint.begin(),matrix_blueprint.end());
    }

    void solve(double f, double gD, double gN, bool fast=false){
        create_blueprint();
        Source(f);
        if (fast){
            fast_boundary(gD,gN);
        }
        else{
            accurate_boundary(gD,gN);
        }
        read_blueprint();

        VectorXd solver_result(N_internal);
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(A);

        if(solver.info()!=Success) {
            cout << "Solver decomposition failed!" << endl;
            return;
        }
        solver_result = solver.solve(rhs);
        if(solver.info()!=Success) {
            cout << "Solver solving failed!" << endl;
            return;
        }
        cout << "Solving succes!" << endl;

        for (int i = 0; i < N_total; i++)
        {
            if (node_grid(i).type==internal_node){
                result(i) = solver_result(node_grid(i).index);
            }
            else{
                result(i) = node_grid(i).eval(gD);
            }
        }    
    } 
// Python conversion
    vector<node_result> get_result(){
        vector<node_result> data(N_total);
        for (int i = 0; i < N_total; i++)
        {
            data[i] = pair(node_grid(i).get_info(),result(i));
        }
        return data;
    }
};



#endif //FDM_HPP