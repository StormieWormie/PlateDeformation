#ifndef MESH_HPP
#define MESH_HPP

#include "setup.hpp"
#include "cluster.hpp"
#include "element.hpp"

class mesh{
public:
//Attributes
    int Number_of_elements;
    vector<element> element_list;
    cluster node_grid;
//Constructors
    mesh(int N_):
            Number_of_elements(N_),
            element_list(N_){}
//Print statements
    void print(){
        cout << "Number of elements: " << Number_of_elements << endl;
        for (element e : element_list){
            e.print();
        }
    }

//Python
    vector<python_node> get_nodes(){
        return node_grid.get_nodes();
    }

    vector<python_element> get_elements(){
        vector<python_element> data(Number_of_elements);
        for (int i = 0; i < Number_of_elements; i++){
            data[i] = element_list[i].get_info();
        }
        return data;
    }
};

class simple_mesh : public mesh{
public:
//Constructors
    simple_mesh(int N):mesh(2*(N-1)*(N-1)){
        node_grid = simple_cluster(N);

        int index(0);
        int a,b,c,d;
        /*
        I look at every 4 nodes in a square grid to create 2 triangular elements
        c---d
        |\  |
        | \ |
        |  \|
        a---b
        */
        for (int i = 0; i < N-1; i++){
            for (int j = 0; j < N-1; j++){
                a = i*N+j;
                b = a+1;
                c = a+N;
                d = c+1;
                element_list[ 2*(i*(N-1)+j) ] = element(node_grid(a),node_grid(b),node_grid(c),index++);
                element_list[2*(i*(N-1)+j)+1] = element(node_grid(d),node_grid(c),node_grid(b),index++);
            }
        }
        /*
        Note, this is a bad meshing as the first and last elements have all nodes on the boundary.
        Instead of adding conditional statements, I simply overwrite the first 2 and last 2 elements by flipping them
        c---d
        |  /|
        | / |
        |/  |
        a---b
        */
        a = 0;b=1;c=N;d=N+1;
        element_list[0] = element(node_grid(c),node_grid(a),node_grid(d),0);
        element_list[1] = element(node_grid(b),node_grid(d),node_grid(a),1);
        d = N*N-1;
        c=d-1;b=d-N;a=b-1;
        index = Number_of_elements-1;
        element_list[index] = element(node_grid(c),node_grid(a),node_grid(d),index);
        index --;
        element_list[index] = element(node_grid(b),node_grid(d),node_grid(a),index);
    }
};

#endif //MESH_HPP