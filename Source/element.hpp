#ifndef ELEMENT_HPP
#define ELEMENT_HPP

#include "setup.hpp"
#include "node.hpp"

enum element_type{boundary_element=1,internal_element=0,untyped_element=-1};
string type2string(element_type input_type){
    switch (input_type)
    {
    case boundary_element:
        return "boundary";
        break;
    case internal_element:
        return "internal";
        break;
    case untyped_element:
        return "untyped";
        break;
    default:
        return "Something went wrong!";
        break;
    }
}

class element{
public:
//Attributes
    int index;                                      //Label this element
    vector<int> node_indices, boundary_indices;     //node labels and registring which nodes are boundary nodes
    int boundary_edge;                              //noting which edge is on the boundary, ab=0,bc=1,ac=2
    double edge_length;                             //The length of the edge, used in boundary integrals
    Matrix2d L,K,KKT;                               //linear transformation, inverse of L and a commonly used product
    Vector2d tau;                                   //seems weird to leave out, is unused but only if f,gD and gN are constants
    double D;                                       //The determinant of L, will come back for every volume integral
    element_type type;                              //type registration to allow easy face value recognition of boundary elements

    element(int index_=-1):
            index(index_),
            node_indices(3),
            boundary_indices(2),
            boundary_edge(-1),
            type(untyped_element)
            {}
    element(node a, node b, node c, int index_=-1):
            index(index_),
            node_indices({a.index,b.index,c.index}),
            boundary_edge(-1){
        tau << a[0], a[1];
        L <<    b[0]-a[0], c[0]-a[0],
                b[1]-a[1], c[1]-a[1];
        D = L(0,0)*L(1,1) - L(0,1)*L(1,0);
        K <<    L(1,1), -L(0,1),
                -L(1,0), L(0,0);
        K /= D;
        KKT = K*K.transpose();
        type = internal_element;
        if(a.type==boundary_node && b.type==boundary_node){ boundary_edge = 0;
                                                            type=boundary_element;
                                                            edge_length = a.distance(b);
                                                            boundary_indices = {a.index, b.index};}
        if(b.type==boundary_node && c.type==boundary_node){ boundary_edge = 1;
                                                            type=boundary_element;
                                                            edge_length = b.distance(c);
                                                            boundary_indices = {b.index, c.index};}
        if(a.type==boundary_node && c.type==boundary_node){ boundary_edge = 2;
                                                            type=boundary_element;
                                                            edge_length = a.distance(c);
                                                            boundary_indices = {a.index, c.index};}
    }   
//Operators
    int operator[](int i){
        return node_indices[i];
    }         
    int operator()(int i){
        return boundary_indices[i];
    }
//Functions
    Matrix<double,10,10> get_hermite_coeffs(){
        Matrix<double,10,10> coeffs;
        coeffs <<   0,0,0,0,27,0,0,-27,-27,0,
                    1,0,0,-3,-13,-3,2,13,13,2,
                    0,L(0,0),L(0,1),-2*L(0,0),-3*(L(0,0)+L(0,1)),-2*L(0,1),L(0,0),3*L(0,0)+2*L(0,1),2*L(0,0)+3*L(0,1),L(0,1),
                    0,L(1,0),L(1,1),-2*L(1,0),-3*(L(1,0)+L(1,1)),-2*L(1,1),L(1,0),3*L(1,0)+2*L(1,1),2*L(1,0)+3*L(1,1),L(1,1),
                    0,0,0,3,-7,0,-2,7,7,0,
                    0,0,0,-L(0,0),2*L(0,0)-L(0,1),0,L(0,0),2*(L(0,1)-L(0,0)),L(0,1)-2*L(0,0),0,
                    0,0,0,-L(1,0),2*L(1,0)-L(1,1),0,L(1,0),2*(L(1,1)-L(1,0)),L(1,1)-2*L(1,0),0,
                    0,0,0,0,-7,3,0,7,7,-2,
                    0,0,0,0,2*L(0,1)-L(0,0),-L(0,1),0,L(0,0)-2*L(0,1),2*(L(0,0)-L(0,1)),L(0,1),
                    0,0,0,0,2*L(1,1)-L(1,0),-L(1,1),0,L(1,0)-2*L(1,1),2*(L(1,0)-L(1,1)),L(1,1);
        return coeffs;
    }
//Print statements
    void print(){
        cout << left
        <<"element: "
        <<setw(10)<< index 
        << "Type: " 
        <<setw(12)<< type2string(type)
        << "nodes: "
        <<right<<setw(8)<<node_indices[0] << ", " 
        <<setw(4) << node_indices[1] << ", " <<setw(4) << node_indices[2] << endl;
    }
//Python
    python_element get_info(){
        python_element data(make_pair(node_indices,make_pair(index,type)));
        return data;
    }
};

#endif //ELEMENT_HPP