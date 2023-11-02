/*
In this document, I will write all setup that can be included as a headerfile.
This includes
    - #include libs
    - typedefs
    - general functions
*/
#ifndef SETUP_HPP
#define SETUP_HPP

//Includes
#include <iostream>
#include <iomanip>
#include <vector>
#include <math.h>
#include <tuple>
#include <utility>

#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

//Namespaces
using namespace std;
using namespace Eigen;

//Typedefs
typedef Eigen::Triplet<double> triplet;
typedef pair<pair<double,double>,pair<int,int>> python_node;
typedef pair<vector<int>,pair<int,int>> python_element;
typedef pair<python_node, double> node_result;
typedef pair<python_node, pair<double, double>> node_result2;
typedef pair<python_node, pair<double, pair<double,double>>> node_result3;

typedef double input_function(double x, double y);

typedef pair<python_node, pair<vector<double>,double>> convergence_result;
//Functions



#endif //SETUP_HPP
