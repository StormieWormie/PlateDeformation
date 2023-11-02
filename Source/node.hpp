#ifndef NODE_HPP
#define NODE_HPP

#include "setup.hpp"

enum node_type{boundary_node=1,internal_node=0,untyped_node=-1};
string type2string(node_type input_type){
    switch (input_type)
    {
    case boundary_node:
        return "boundary";
        break;
    case internal_node:
        return "internal";
        break;
    case untyped_node:
        return "untyped";
        break;
    default:
        return "Something went wrong!";
        break;
    }
}

struct node
{
public:
//Atributes
    int index;
    node_type type;
    Vector2d position;
//Constructors
    node(Vector2d position_, int index_=-1, node_type type_=untyped_node):
        position(position_),
        index(index_),
        type(type_){}
    node(int index_=-1, node_type type_=untyped_node):
        index(index_),
        type(type_){}
//Operators
    double operator[](int i){
        return position(i);
    }
//Functions
    double distance(node other){
        Vector2d d(other.position - position);
        return d.norm();
    }
    double eval(input_function func){
        return func(position(0),position(1));
    }
    double eval(double val){
        return val;
    }
//Routines
    void set_index(int index_){index = index_;}
    void set_type(node_type type_){type = type_;}
    void set_position(Vector2d position_){position = position_;}
    void set_type_and_index(node_type type_, int index_){set_type(type_);set_index(index_);}
    
//Print statements
    void print(){
        cout << left
        <<"node: "
        <<setw(10)<< index 
        << "Type: " 
        <<setw(12)<< type2string(type)
        << "Position: "
        <<right<<setw(8)<<position(0) << ", " <<left
        <<setw(8)<< position(1) << endl;
    }

//Python
    python_node get_info(){
        python_node data(make_pair(position(0),position(1)),make_pair(index,type));
        return data;
    }
};


#endif //NODE_HPP