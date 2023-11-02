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
    - 2nd order accurate (Neumann condition approximated using central difference)
