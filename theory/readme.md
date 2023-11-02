# Theory
## The plate equation
As I want to be able to solve 4th order PDE's, the first testcase is one such PDE with no additional terms to simply get to the point.
$$\nabla^4u(x,y) = f(x,y) \qquad \text{for } (x,y)\in \Omega$$
For boundary condition, I choose clamped boundary conditions. A combination of both Dirichlet and Neumann conditions.
$$u(x,y) = g_1(x,y) \qquad \text{for } (x,y)\in \Gamma$$
$$\frac{\partial u}{\partial \vec{n}}|_{x,y} = g_2(x,y) \qquad \text{for } (x,y)\in \Gamma$$

## The weak form
Using testfunction $v$, I write the weak form of the system.
