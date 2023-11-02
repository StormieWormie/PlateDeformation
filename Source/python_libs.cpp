#ifndef PYTHON
#define PYTHON

#include "setup.hpp"
#include "FiniteDifferenceMethod.hpp"
#include "FiniteElementMethod.hpp"
#include "Convergence.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

namespace py = pybind11;

void foo(){
    cout << "The empty function can be called." << endl;
}


PYBIND11_MODULE(FiniteDifferenceMethod, handle){
    handle.doc() = "Finite Difference Method";
    handle.def("foo", foo, "Just an empty function to be called");

    py::class_<simple_cluster>(handle, "cluster")
        .def(py::init<int>())
        .def("print",&simple_cluster::print)
        .def("get_nodes",&cluster::get_nodes);
    
    py::class_<FiniteDifferenceMethod>(handle, "FiniteDifferenceMethod")
        .def(py::init<int>())
        .def("solve",&FiniteDifferenceMethod::solve)
        .def("get_result",&FiniteDifferenceMethod::get_result);
}

PYBIND11_MODULE(FiniteElementMethod, handle){
    handle.doc() = "Finite Element Method";
    handle.def("foo", foo, "Just an empty function to be called");

    py::class_<simple_mesh>(handle, "mesh")
        .def(py::init<int>())
        .def("print",&simple_mesh::print)
        .def("get_nodes",&simple_mesh::get_nodes)
        .def("get_elements",&simple_mesh::get_elements);

    py::class_<MixedFiniteElementMethod>(handle, "MixedFiniteElementMethod")
        .def(py::init<int>())
        .def("solve",&MixedFiniteElementMethod::solve)
        .def("get_result",&MixedFiniteElementMethod::get_result);
    
    py::class_<ExoticFiniteElementMethod>(handle, "ExoticFiniteElementMethod")
        .def(py::init<int>())
        .def("solve",&ExoticFiniteElementMethod::solve)
        .def("get_result",&ExoticFiniteElementMethod::get_result);
}

PYBIND11_MODULE(Convergence, handle){
    handle.doc() = "Convergence Analysis";

    py::class_<Convergence>(handle, "ConvergenceAnalysis")
        .def(py::init<int>())
        .def(py::init<int,int,int>())
        .def(py::init<int,int,int,int>())
        .def("get_N",&Convergence::get_N)
        .def("compute_fdm",&Convergence::compute_fdm)
        .def("compute_fem",&Convergence::compute_mixed_fem)
        .def("compute_norms",&Convergence::compute_norms)
        .def("get_data",&Convergence::get_convergence_data)
        .def("get_norms",&Convergence::get_norms)
        ;

    py::class_<RConvergence>(handle, "RConvergenceAnalysis")
        .def(py::init<int>())
        .def(py::init<int,int>())
        .def(py::init<int,int,int>())
        .def("get_N",&RConvergence::get_N)
        .def("compute_fdm",&RConvergence::compute_fdm)
        .def("get_data",&RConvergence::get_convergence_data)
        ;
}
#endif //PYTHON