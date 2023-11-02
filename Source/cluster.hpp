#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include "setup.hpp"
#include "node.hpp"

class cluster
{
public:
//Attributes
    int Number_of_nodes;
    vector<node> node_list;
//Constructors
    cluster(int Number_of_nodes_=0):
        Number_of_nodes(Number_of_nodes_),
        node_list(Number_of_nodes_){}
//Operator
    node operator()(int i){
        return node_list[i];
    }
//Print statements
    void print(){
        cout << "Number of nodes: " << Number_of_nodes << endl;
        for (node n : node_list){
            n.print();
        }   
    }

//Python conversion
    vector<python_node> get_nodes(){
        vector<python_node> data(Number_of_nodes);
        for (int i = 0; i < Number_of_nodes; i++)
        {
            data[i] = node_list[i].get_info();
        }
        return data;
    }
};

class simple_cluster : public cluster
{
public:
//Constructors
    simple_cluster(int N):cluster(N*N){
        double h = 1./double(N-1);
        int index(0);
        for (int j = 0; j < N; j++){
            for (int i = 0; i < N; i++){
                node_list[N*j+i].set_position({i*h,j*h});
            }
        }
        for (int j = 1; j < N-1; j++){
            for (int i = 1; i < N-1; i++){
                node_list[N*j+i].set_type_and_index(internal_node, index++);
            }
        }
        for (int i = 0; i < N-1; i++){node_list[i].set_type_and_index(boundary_node, index++);}
        for (int i = 0; i < N-1; i++){node_list[N-1+i*N].set_type_and_index(boundary_node, index++);}
        for (int i = 0; i < N-1; i++){node_list[Number_of_nodes-1-i].set_type_and_index(boundary_node, index++);}
        for (int i = 0; i < N-1; i++){node_list[Number_of_nodes-N-N*i].set_type_and_index(boundary_node, index++);}
        

    }
};

#endif //CLUSTER_HPP