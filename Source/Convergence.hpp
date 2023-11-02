#ifndef CONVERGENCE_HPP
#define CONVERGENCE_HPP

#include "setup.hpp"
#include "FiniteDifferenceMethod.hpp"
#include "FiniteElementMethod.hpp"

class Convergence
{
public:
//Attributes
    int Nstart, factor, hyperstep, Nsmall;
    vector<int> N;
    int Nhyper, hyperfactor;

    vector<convergence_result> convergence_data;
    vector<double> L2Norm, MaxNorm;
//Constructors
    Convergence(int Nstart_, int increase_=2, int hyperstep_= 2, int Nsmall_=3):
            Nstart(Nstart_),
            factor(increase_),
            hyperstep(hyperstep_),
            Nsmall(Nsmall_),
            N(Nsmall)
        {
        N[0] = Nstart;
        hyperfactor = 1;
        for (int i = 1; i < Nsmall; i++)
        {
            N[i] = factor * (N[i-1]-1) + 1;
            hyperfactor *= factor;
        }
        Nhyper = N[Nsmall-1];
        for (int i = 0; i < hyperstep; i++)
        {
            Nhyper = factor*(Nhyper-1)+1;
            hyperfactor *= factor;
        }
    }

    vector<int> get_N(){
        vector<int> data(Nsmall+1);
        for (int i = 0; i < Nsmall; i++)
        {
            data[i] = N[i];
        }
        data[Nsmall] = Nhyper;
        return data;
    }

    void compute_fdm(double f, double Dirichlet, double Neumann, bool fast = false){
        vector<FiniteDifferenceMethod> drivers(Nsmall+1);
        for (int i = 0; i < Nsmall; i++)
        {
            drivers[i] = FiniteDifferenceMethod(N[i]);
        }
        
        drivers[Nsmall] = FiniteDifferenceMethod(Nhyper);
        
        for (int i = 0; i < Nsmall+1; i++)
        {
            drivers[i].solve(f,Dirichlet,Neumann,fast);
        }

        python_node node_info;
        vector<double> results(Nsmall);
        double hyper_result;
        convergence_data.resize(drivers[0].N_total);
        int scale;

        for (int i = 0; i < N[0]; i++){
            for (int j = 0; j < N[0]; j++){
                node_info = drivers[0].node_grid(i*N[0]+j).get_info();
                scale = 1;
                for (int k = 0; k < Nsmall; k++)
                {
                    results[k] = drivers[k].result(scale*(i*N[k]+j));
                    scale *= factor;
                }
                
                hyper_result = drivers[Nsmall].result(hyperfactor*(i*Nhyper+j));
                convergence_data[i*N[0]+j] = pair(node_info,pair(results,hyper_result));
            }
        }        
    }

    void compute_mixed_fem(double f, double Dirichlet, double Neumann){
        vector<MixedFiniteElementMethod> drivers(Nsmall+1);
        for (int i = 0; i < Nsmall; i++)
        {
            drivers[i] = MixedFiniteElementMethod(N[i]);
        }
        drivers[Nsmall] = MixedFiniteElementMethod(Nhyper);
        for (int i = 0; i < Nsmall+1; i++)
        {
            drivers[i].solve(f,Dirichlet,Neumann);
        }

        python_node node_info;
        vector<double> results(Nsmall);
        double hyper_result;
        convergence_data.resize(drivers[0].N_total);
        int scale, N_total;

        for (int i = 0; i < N[0]; i++){
            for (int j = 0; j < N[0]; j++){
                node_info = drivers[0].element_grid.node_grid(i*N[0]+j).get_info();
                scale = 1;
                for (int k = 0; k < Nsmall; k++)
                {
                    N_total = drivers[k].N_total;
                    results[k] = drivers[k].result(N_total + scale*(i*N[k]+j));
                    scale *= factor;
                }
                
                N_total = drivers[Nsmall].N_total;
                hyper_result = drivers[Nsmall].result(N_total + hyperfactor*(i*Nhyper+j));
                convergence_data[i*N[0]+j] = pair(node_info,pair(results,hyper_result));
            }
        }
        
        
    }

    void compute_norms(){
        double error;
        MaxNorm = {0,0,0};
        L2Norm = {0,0,0};
        for (convergence_result result : convergence_data){            
            for (int i = 0; i < 3; i++){
                error = result.second.first[i] - result.second.second;
                error = abs(error);
                L2Norm[i] += error*error;
                if(error > MaxNorm[i]){
                    MaxNorm[i] = error;
                }
            }
        }
        for (int i = 0; i < 3; i++)
        {
            L2Norm[i] = sqrt(L2Norm[i]);
        }
    }

    vector<convergence_result> get_convergence_data(){
        return convergence_data;
    }
    pair<vector<double>,vector<double>> get_norms(){
        return pair(L2Norm,MaxNorm);
    }
};

typedef pair<python_node, vector<double>> Rconvergence_result;

class RConvergence
{
public:
//Attributes
    int Nstart, factor, Nsmall;
    vector<int> N;

    vector<Rconvergence_result> convergence_data;
    vector<double> L2Norm, MaxNorm;
//Constructors
    RConvergence(int Nstart_ = 5, int factor_ = 2, int Nsmall_ = 3):
            Nstart(Nstart_),factor(factor_),Nsmall(Nsmall_),
            N(Nsmall),L2Norm(Nsmall),MaxNorm(Nsmall)
    {
        N[0] = Nstart;
        for (int i = 1; i < Nsmall; i++)
        {
            N[i] = factor*(N[i-1]-1) + 1;
        }
    }

    vector<int> get_N(){return N;}

    void compute_fdm(double Source, double Dirichlet, double Neumann, bool fast = false){
        vector<FiniteDifferenceMethod> drivers(Nsmall);
        for (int i = 0; i < Nsmall; i++)
        {
            drivers[i] = FiniteDifferenceMethod(N[i]);
            drivers[i].solve(Source,Dirichlet,Neumann,fast);
        }
        python_node node_info;
        vector<double> results(Nsmall);
        convergence_data.resize(drivers[0].N_total);

        int scale;
        for (int i = 0; i < N[0]; i++){
            for (int j = 0; j < N[0]; j++){
                node_info = drivers[0].node_grid(i*N[0]+j).get_info();
                scale = 1;
                for (int k = 0; k < Nsmall; k++)
                {
                    results[k] = drivers[k].result(scale*(i*N[k]+j));
                    scale *= factor;
                }
                convergence_data[i*N[0]+j] = pair(node_info, results);
            }
        }
    }



    vector<Rconvergence_result> get_convergence_data(){
        return convergence_data;
    }
};


#endif //CONVERGENCE_HPP   