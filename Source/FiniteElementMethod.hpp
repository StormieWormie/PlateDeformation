#ifndef FEM_HPP
#define FEM_HPP

#include "setup.hpp"
#include "mesh.hpp"


class FiniteElementMethod{
public:
//Attributes
    simple_mesh element_grid;
    int N_total, N_internal, N_boundary;
    int Domain_size, Problem_size;
    vector<triplet> mass_blueprint, stress_blueprint;
    SparseMatrix<double> mass_matrix, stress_matrix;


    vector<triplet> matrix_blueprint;
    SparseMatrix<double> A; VectorXd rhs;
    VectorXd result;
//Constructors
    FiniteElementMethod(int N):
        element_grid(N),
        N_total(N*N),N_internal((N-2)*(N-2)),N_boundary(N_total-N_internal){}

//Routines
    void resize(){
        A.resize(Domain_size,Domain_size);
        rhs.resize(Domain_size);
        result.resize(Domain_size);
        rhs.setZero();
        result.setZero();
    }
    void read_blueprint(){
        A.setFromTriplets(matrix_blueprint.begin(),matrix_blueprint.end());
    }
    void read_mass_blueprint(){mass_matrix.setFromTriplets(mass_blueprint.begin(),mass_blueprint.end());}
    void read_stress_blueprint(){stress_matrix.setFromTriplets(stress_blueprint.begin(),stress_blueprint.end());}
    
};

class MixedFiniteElementMethod : public FiniteElementMethod{
public:
//Attributes

//Constructors
    MixedFiniteElementMethod(int N=5):
            FiniteElementMethod(N)
            {
        Domain_size = 2*N_total;                //you solve for 2 functions
        Problem_size = N_total + N_internal;    //One function has Dirichlet conditions and needs less solving
        resize();
        mass_matrix.resize(N_total,N_total);
        stress_matrix.resize(N_total,N_total);
    }

    void create_mass_blueprint(){
        Matrix3d Coeffs, Integrals, reference_matrix, element_matrix;
        int ii,jj;
        Coeffs <<   1,-1,-1,
                    0, 1, 0,
                    0, 0, 1;
        Integrals <<    1./2.,1./6.,1./6.,
                        1./6.,1./12.,1./24.,
                        1./6.,1./24.,1./12.;
        reference_matrix = Coeffs*Integrals*Coeffs.transpose();
        for (element e : element_grid.element_list){
            element_matrix = e.D * reference_matrix;
            for (int i = 0; i < 3; i++){
                ii = e[i];
                for (int j = 0; j < 3; j++){
                    jj = e[j];
                    mass_blueprint.push_back(triplet(ii,jj,element_matrix(i,j)));
                }
            }
        }
    }
    void create_stress_blueprint(){
        Matrix3d Coeffs, element_matrix;
        Matrix<double, 2, 3> gradients;
        int ii,jj;
        Coeffs <<   1,-1,-1,
                    0, 1, 0,
                    0, 0, 1;
        gradients <<    -1,1,0,
                        -1,0,1;
        
        for (element e : element_grid.element_list){
            element_matrix = .5* gradients.transpose() * e.KKT * gradients * e.D;
            for (int i = 0; i < 3; i++){
                ii = e[i];
                for (int j = 0; j < 3; j++){
                    jj = e[j];
                    stress_blueprint.push_back(triplet(ii,jj,element_matrix(i,j)));
                }
            }
        }
    }
    void create_blueprint(){
        for (triplet T : mass_blueprint){
            matrix_blueprint.push_back(triplet(T.row(),T.col(),-T.value()));
        }
        for (triplet T : stress_blueprint){
            matrix_blueprint.push_back(triplet(T.row(),T.col()+N_total,T.value()));
            matrix_blueprint.push_back(triplet(T.col()+N_total,T.row(),T.value()));
        }
    }

    void add_neumann(double Neumann){
        Vector2d element_vector, reference_vector;
        reference_vector << .5,.5;
        reference_vector *= Neumann;
        for (element e : element_grid.element_list){
            if (e.type==boundary_element){
                element_vector = reference_vector*e.edge_length;
                for (int i = 0; i < 2; i++){
                    rhs(e(i)) += element_vector(i);
                }
            }
        }
    }

    void solve(double f, double Dirichlet=0, double Neumann=0){
        create_mass_blueprint();
        read_mass_blueprint();

        create_stress_blueprint();
        read_stress_blueprint();

        create_blueprint();
        read_blueprint();

        VectorXd f_vec(N_total);
        f_vec.setConstant(f);
        rhs.segment(N_total,N_internal) += mass_matrix.topRows(N_internal) * f_vec;

        VectorXd ub(N_boundary);
        ub.setConstant(Dirichlet);
        rhs.head(Problem_size) -= A.topRightCorner(Problem_size,N_boundary) * ub;

        add_neumann(Neumann);

        SparseMatrix<double> problem_matrix(Problem_size,Problem_size);
        problem_matrix = A.topLeftCorner(Problem_size,Problem_size);

        VectorXd problem_vector(Problem_size);
        problem_vector = rhs.head(Problem_size);

        VectorXd solver_result(N_internal);
        SparseLU<SparseMatrix<double>> solver;
        solver.compute(problem_matrix);

        if(solver.info()!=Success) {
            cout << "Solver decomposition failed!" << endl;
            return;
        }
        solver_result = solver.solve(problem_vector);
        if(solver.info()!=Success) {
            cout << "Solver solving failed!" << endl;
            return;
        }
        cout << "Solving succes!" << endl;

        for (int i = 0; i < N_total; i++){
            result(i) = solver_result(element_grid.node_grid(i).index);
            if(element_grid.node_grid(i).type==internal_node){
                result(N_total+i) = solver_result(N_total + element_grid.node_grid(i).index);
            }
            else{
                result(N_total+i) = Dirichlet;
            }
        }
    }

//Python
    vector<node_result2> get_result(){
        vector<node_result2> data(N_total);
        for (int i = 0; i < N_total; i++){
            data[i] = make_pair(element_grid.node_grid(i).get_info(),make_pair(result(i),result(N_total+i)));
        }
        return data;
    }
};


class ExoticFiniteElementMethod : public FiniteElementMethod{
public:
//Constructors
    ExoticFiniteElementMethod(int N):
            FiniteElementMethod(N)
            {
        Domain_size = 3*N_total + element_grid.Number_of_elements;                //you solve for 3 functions and add nodes in the middle of elements
        Problem_size = 3*N_internal + element_grid.Number_of_elements;            //all functions have Dirichlet conditions
        resize();
        mass_matrix.resize(Domain_size,Domain_size);
    }

//Functions
    vector<int> get_indices(element e){
        vector<int> indices(10);
        indices[0] = e.index;
        for (int i = 0; i < 3; i++){
            for (int j = 0; j < 3; j++){
                indices[1+3*i+j] = element_grid.Number_of_elements + 3*e[i] + j;
            }
        }
        return indices;
    }

//Routines
    void create_mass_blueprint(){
        Matrix<double,10,10> Coeffs, Integrals, element_matrix;
        int ii,jj;
        vector<int> index_list;
        Integrals << 
            2,6,6,12,24,12,20,60,60,20,
            6,12,24,20,60,60,30,120,180,120,
            6,24,12,60,60,20,120,180,120,30,
            12,20,60,30,120,180,42,210,420,420,
            24,60,60,120,180,120,210,420,420,210,
            12,60,20,180,120,30,420,420,210,42,
            20,30,120,42,210,420,56,336,840,1120,
            60,120,180,210,420,420,336,840,1120,840,
            60,180,120,420,420,210,840,1120,840,336,
            20,120,30,420,210,42,1120,840,336,56;
        for (int i = 0; i < 10; i++){
            for (int j = 0; j < 10; j++){
                    Integrals(i,j) = 1./Integrals(i,j);
                }
            }
        for (element e : element_grid.element_list){
            Coeffs = e.get_hermite_coeffs();
            element_matrix = Coeffs*Integrals*Coeffs.transpose()*e.D;
            index_list = get_indices(e);
            for (int i = 0; i < 10; i++){
                ii = index_list[i];
                for (int j = 0; j < 10; j++){
                    jj = index_list[j];
                    mass_blueprint.push_back(triplet(ii,jj,element_matrix(i,j)));
                }
            }
        }
    }

    void create_matrix_blueprint(){
        Matrix<double,10,10> Coeffs, Integrals, reference_matrix, element_matrix;
        Matrix<double,10,1> a0,ax,ay;
        Matrix<double,10,10> b00,b0x,b0y,bxx,bxy,byy;
        int ii,jj;
        vector<int> index_list(10);
        for (element e : element_grid.element_list){
            a0 << 0,0,0,2*e.KKT(0,0),e.KKT(1,0)+e.KKT(0,1),2*e.KKT(1,1),0,0,0,0;
            ax << 0,0,0,0,0,0,3*a0(3),2*a0(4),a0(5),0;
            ay << 0,0,0,0,0,0,0,a0(3),2*a0(4),3*a0(5);
            b00 = a0*a0.transpose();
            b0x = a0*ax.transpose() + ax*a0.transpose();
            b0y = a0*ay.transpose() + ay*a0.transpose();
            bxx = ax*ax.transpose();
            bxy = ax*ay.transpose() + ay*ax.transpose();
            byy = ay*ay.transpose();
            Integrals = 1./2.*b00 +
                        1./6.*(b0x+b0y) +
                        1./12.*(bxx + .5*bxy + byy);
            Coeffs = e.get_hermite_coeffs();
            element_matrix =Coeffs * Integrals * Coeffs.transpose() * e.D;
            index_list = get_indices(e);
            for (int i = 0; i < 10; i++){
                ii = index_list[i];
                for (int j = 0; j < 10; j++){
                    jj = index_list[j];
                    matrix_blueprint.push_back(triplet(ii,jj,element_matrix(i,j)));
                }
            }
        }
    }

    void solve(double f, double Dirichlet=0, double Neumann=0){
        create_mass_blueprint();
        read_mass_blueprint();

        create_matrix_blueprint();
        read_blueprint();

        VectorXd f_vec(Domain_size);
        f_vec.setZero();
        for (int i = 0; i < element_grid.Number_of_elements; i++){
            f_vec(i) = f;
        }
        for (int i = 0; i < N_total; i++){
            f_vec(element_grid.Number_of_elements+3*i) = f;
        }
        rhs = mass_matrix*f_vec;
        
        VectorXd problem_vector(Problem_size), solver_result(Problem_size);
        problem_vector = rhs.head(Problem_size);

        
        VectorXd boundary_vector(Domain_size);
        boundary_vector.setZero();
        
        for (node n : element_grid.node_grid.node_list){
            if(n.type == boundary_node){
                boundary_vector(element_grid.Number_of_elements+3*n.index) = Dirichlet;
                if(n[0]==1){boundary_vector(element_grid.Number_of_elements+3*n.index+1) = Neumann;}
                if(n[0]==0){boundary_vector(element_grid.Number_of_elements+3*n.index+1) = -Neumann;}
                
                if(n[1]==1){boundary_vector(element_grid.Number_of_elements+3*n.index+2) = Neumann;}
                if(n[1]==0){boundary_vector(element_grid.Number_of_elements+3*n.index+2) = -Neumann;}
                
            }
        }
        problem_vector -= A.topRows(Problem_size)*boundary_vector;
        

        SparseMatrix<double> problem_matrix(Problem_size,Problem_size);
        problem_matrix = A.topLeftCorner(Problem_size,Problem_size);
        

        SparseLU<SparseMatrix<double>> solver;
        solver.compute(problem_matrix);

        if(solver.info()!=Success) {
            cout << "Solver decomposition failed!" << endl;
            return;
        }
        solver_result = solver.solve(problem_vector);
        if(solver.info()!=Success) {
            cout << "Solver solving failed!" << endl;
            return;
        }
        cout << "Solving succes!" << endl;

        result.head(element_grid.Number_of_elements) = solver_result.head(element_grid.Number_of_elements);
        for (int i = 0; i < N_total; i++){
            if(element_grid.node_grid(i).type==internal_node){
                for (int k = 0; k < 3; k++){
                    result(element_grid.Number_of_elements + 3*i + k) = solver_result(element_grid.Number_of_elements + 3*element_grid.node_grid(i).index + k) + 
                                                                        boundary_vector(element_grid.Number_of_elements + 3*element_grid.node_grid(i).index + k);
                }
            }
            else{
                for (int k = 0; k < 3; k++){
                    result(element_grid.Number_of_elements + 3*i + k) = boundary_vector(element_grid.Number_of_elements + 3*element_grid.node_grid(i).index + k);
                }
                
                //result(element_grid.Number_of_elements + 3*i) = boundary_vector(element_grid.Number_of_elements + 3*element_grid.node_grid(i).index);
                //result(element_grid.Number_of_elements + 3*i + 1) = 0; //incorrect
                //result(element_grid.Number_of_elements + 3*i + 2) = 0; //incorrect
            }
        }
    }

    

    vector<node_result3> get_result(){
        vector<node_result3> data(N_total);
        int I;
        for (int i = 0; i < N_total; i++)
        {
            I = element_grid.Number_of_elements + 3*i;
            data[i] = make_pair(element_grid.node_grid(i).get_info(), make_pair(result(I), make_pair(result(I+1),result(I+2))));
        }
        return data;
        
    }
};

#endif //FEM_HPP