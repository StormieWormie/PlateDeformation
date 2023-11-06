# Solving fourth order differential equations
In this repository, I build the tools to solve fourth order differential equations using the finite element method.

## An introduction
For my bachelor thesis project, I attempted to implement the finite element method to simulate self-induced curvature of biological membranes in C++. This system was modeled using a coupled PDE system consisting of two energy sources that needed to be minimized. The assumption is that freely moving curved proteins on/in the membrane cause the curvature of membranes. 
- the Canham-Helfrich equation (Surface tension)
- the Cahn-Hilliard equation (Phase separation of the curvature inducing proteins)

Both equations are challenging in their own right but my interest primarily stems from the former component. The surface tension is described with a 4th order differential equation, something that I do not advice to attempt to solve for your first finite element method implementation. At the time, however, I did not realize this and this project was my first contact with FEM.
Depending on the properties of the curvature inducing proteins, the latter equation can be non-linear. For my thesis, the property of the proteins was assumed to be linear and the non-linear component was ignored.

![hello](https://upload.wikimedia.org/wikipedia/commons/thumb/9/9e/Blausen_0350_EndoplasmicReticulum.png/330px-Blausen_0350_EndoplasmicReticulum.png)

Image Source: [Blausen.com staff (2014). "Medical gallery of Blausen Medical 2014". WikiJournal of Medicine 1 (2). DOI:10.15347/wjm/2014.010. ISSN 2002-4436. - Own work, CC BY 3.0](https://en.wikipedia.org/wiki/Endoplasmic_reticulum#/media/File:Blausen_0350_EndoplasmicReticulum.png)

Afterwards, in my master programme, I learned how the finite element method is implemented. First, I learned how the finite element method is implemented for second order differential equations. Next, I learned a method to solve the beam equation, a fourth order ODE. To round of this set of courses, I chose to work on an individual project combining the knowledge from different courses to properly solve a fourth order PDE, using object-oriented C++ as a framework. 

## The state of this repo
In short, the state is **ONGOING**
### Currently done
#### Grid generation
- An equidistant grid can be generated with nodes, clearly seperated by an internal and boundary section.
#### Numerical solvers
- The finite difference method can be used with a 13-point stencil to solve the clamped plate equation with constant heterogeneous boundary conditions and source.
    - 1st order accurate (Neumann condition approximated using forward difference)
    - 2nd order accurate (Neumann condition approximated using central difference
- The mixed finite element method can be used to solve the clamped plate equation with linear interpolation functions as a basis
    - 2nd order accurate
#### Quality of life
- Using pybind, the results of the solver classes can be easily read in python, allowing for easier testing and visualization.
### Currently in the works
- The finite element method using cubic hermite interpolation functions as a basis
    - The mass matrix seems to function properly.
    - The solution is simply not consistent with the mixed fem and fdm solutions
        - The mass matrix seems to function well through testing with simple functions, more complex functions testing is possible but easier if functionality of the nodes and grid improves
        - I suspect boundary conditions are not implemented correctly but a more abstract implementation of the nodes/grid eases the testing.
### Future improvements
- Increase node functionality (**now**)
    - Add function class to give to nodes to evaluate those functions and their derivatives on the nodes.
    - Add a more extensive labeling to distinguish not simpley internal/boundary nodes but also have different boundary labels.
- Increase more general grid generation methods
    - Add a method to generate a grid, given a user-input boundary (will still be square at first)
    - Add a simple interface in python to "draw" domains
- Add different dimension elements
    - allowing N-dimensional boundary elements to just use (N-1)-dimensional internal elements
- Create some function mapping
    - Function classes should be definable in python to redduce recompilation.
# Dependencies
- Runtime
    - Eigen3 (https://eigen.tuxfamily.org/index.php?title=Main_Page)
        - Solving the sparse matrix-vector problem at the end
    - Pybind (https://pybind11.readthedocs.io/en/stable/index.html)
        - Linking C++ and python to use the computational power of C++ while using the visualization of python
- Preperatory
    - Numpy
        - What numerical code doesn't use numpy?
    - Sympy (https://www.sympy.org/en/index.html)
        - The FEM uses isoparametric and a reference element, integration and derivation is done beforehand on the reference element and all elements use those solutions as a basis to avoid having to do calculus for each element. Sympy is used to ease the workload (integration basis of cubic polynomial requires integration of 100 functions) and solve abstract problems (6/10 cubic polynomials depend on the position of the element and would require solving a linear problem, using an abstract element and solving for it allows me to just fill in the solution for each element instead of solving it each time.)

# How to use the code
The Source directory contains a CMake file that can be used to make all binaries. The most useful components are
- Make exe (This makes the main.cpp file to "Experiment")
- Make libs (This builds both FiniteDifferenceMethod and FiniteElementMethod, 2 importable modules for python)

# Results
## Convergence
![conv1](https://github.com/StormieWormie/PlateDeformation/assets/46678214/ea5255d6-4030-4e0f-838b-50d485f2b2c6)

## Interpolation
Testing the interpolation $e^{x^2+y^2}$
![interpolation_error](https://github.com/StormieWormie/PlateDeformation/assets/46678214/6e8402e1-2400-4372-b7e4-093d5f062a6d)

Cubic hermite interpolation gives a fourth order accurate approximation. It seems worth investigating further.

# Theory
For different parts of the theory, you can look at the theory directory which looks at different components of the problem.






