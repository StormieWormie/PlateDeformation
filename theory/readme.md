# Theory
## The plate equation
As I want to be able to solve 4th order PDE's, the first testcase is one such PDE with no additional terms to simply get to the point.
$$\nabla^4u(x,y) = f(x,y) \qquad \text{for } (x,y)\in \Omega$$
For boundary condition, I choose clamped boundary conditions. A combination of both Dirichlet and Neumann conditions.
$$u(x,y) = g_1(x,y) \qquad \text{for } (x,y)\in \Gamma$$
$$\frac{\partial u}{\partial \vec{n}}|_{x,y} = g_2(x,y) \qquad \text{for } (x,y)\in \Gamma$$

## The weak form
Using testfunction $v$, I write the weak form of the system.
$$\int_\Omega \nabla^2v\nabla^2u d\vec{x} = \int_\Omega vf d\vec{x} -
\oint_\Gamma v (\nabla^3u\cdot\vec{n}) d\vec{x}
+\oint_\Gamma \nabla^2u(\nabla v\cdot\vec{n})d\vec{x}$$

For the simplest implementations of the finite element method, there is an issue here. Linear interpolation function don't work well with the $\nabla^2$ operator as it will always evaluate to 0. Eliminating terms is nice but only to an extend, here it goes too far and the problem becomes trivial and nonsensical. I look at 2 methods around this issue.
- Mixed finite element method
  - The equation is too hard so I change the equation
- cubic hermite basis functions
  - Linear basis functions don't work but cubic functions might
 
The boundary integrals pose a challenge but if the testfunction has homogeneous clamped plate conditions, both evaluate to 0. 

## Mixed finite element method
The equation can be rewritten as a coupled system.
$$-\nabla^2u(x,y) = w(x,y) \qquad \text{for } (x,y)\in \Omega$$
$$-\nabla^2w(x,y) = f(x,y) \qquad \text{for } (x,y)\in \Omega$$

Rewritten in the weak form.

$$-\int_\Omega vwd\vec{x} + \int_\Omega\nabla v\cdot\nabla ud\vec{x}
    =\oint_\Gamma v  \frac{\partial u}{\partial \vec{n}} d\vec{x} $$
    
$$ \int_\Omega \nabla v_0\cdot\nabla w d\vec{x} =
    \oint_\Gamma v_0 \frac{\partial w}{\partial \vec{n}} d\vec{x} + \int_\Omega v_0fd\vec{x}$$

For symmetry of the matrices, $v_0$ has homogeneous Dirichlet conditions but $v$ does not. Linear interpolation functions can be used for these equations but I cannot eliminate both boundary integrals. 

## Isoparametric elements
For the elements, I use isoparametric elements. This means that all elements are defined as a transformation with respect to a reference element. The functions on each element are also defined as a transformation of the functions on the reference element. This allows me to only have to define the reference functions. 

![transformation](https://github.com/StormieWormie/PlateDeformation/assets/46678214/d5adc4e0-c8a7-49c8-8b31-5d9faae0dfe4)

The transformation from the reference coordinates $ \vec{\xi} $ to the real coordinates can be written as
$$\vec{x} = L\vec{\xi} + \tau$$
