# Mixed finite element method
Here, I go over some of the mathematical preparation for implementing the mixed finite element method
## The basis function
The basis for my mixed fem consists of linear interpolation functions (hatfunctions). On the reference element, the 3 functions can already be defined. The basis of linear functions $\phi_g(\vec{\xi})$ is used here.

$$\phi_g(\vec{\xi}) = \begin{bmatrix}
1\\ 
\xi \\
\eta \end{bmatrix}$$

$$C_1 = \begin{bmatrix}
1&-1&-1\\
0&1&0\\
0&0&1
\end{bmatrix}$$

$$
\phi (\vec{\xi}) =
\begin{bmatrix}
\phi_0(\vec{\xi})\\
\phi_1(\vec{\xi})\\
\phi_2(\vec{\xi})
\end{bmatrix}=
C_1 \phi_g(\vec{\xi})
$$

## The mass matrix
The element mass matrix can be found by computing 
$$\int_{\Omega_e} \phi (\vec{\xi}) \phi (\vec{\xi})^T |L|d\vec{\xi}$$
Here, the reason for the previous components becomes a bit more clear as I can take $C_1$ out of the integral.
$$|L| C_1 \int_{\Omega_e} \phi_g (\vec{\xi}) \phi_g (\vec{\xi})^T d\vec{\xi}C_1^T$$
But the basis of linear functions is already known so the integral can be computed beforehand.

$$|L| C_1 
\begin{bmatrix}
        \frac{1}{2} & \frac{1}{6} & \frac{1}{6} \\
        \frac{1}{6} & \frac{1}{24} & \frac{1}{12} \\
        \frac{1}{6} & \frac{1}{12} & \frac{1}{24}
    \end{bmatrix}
C_1^T$$

And with the linear functions $C_1$ is constant (in contrast to cubic hermite functions) so the entire matrix is even simpler.

$$
|L|
\begin{bmatrix}
        \frac{1}{12} & \frac{1}{24} & \frac{1}{24}\\ 
        \frac{1}{24} & \frac{1}{12} & \frac{1}{24}\\ 
        \frac{1}{24} & \frac{1}{24} & \frac{1}{12}
    \end{bmatrix}
$$

For every element, the only difference in the element mass matrix is a scalar $|L|$. 
## The stiffness matrix
The element stiffness matrix can be found by computing 
$$\int_{\Omega_e} \nabla_{\mathbf{x}}\phi_i (\vec{\xi}) \cdot \nabla_{\mathbf{x}}\phi_j (\vec{\xi}) |L|d\vec{\xi}$$

$$ \int_{\Omega_e} (\nabla_{\mathbf{\xi}} \psi_i)^T K K^T\nabla_{\mathbf{\xi}} \psi_j |L|d\xi d\eta$$

The gradient of the linear functions can be computed beforehand and filled in.

$$\int_{\Omega_e}\begin{bmatrix}
-1&-1\\
1&0\\
0&1
\end{bmatrix}
KK^T
\begin{bmatrix}
-1&1&0\\
-1&0&1
\end{bmatrix}
|L|d\xi d\eta$$

There are multiple linear operations but no more variables. I can integrate beforehand.

$$
\frac{1}{2}
|L|
\begin{bmatrix}
-1&-1\\
1&0\\
0&1
\end{bmatrix}
KK^T
\begin{bmatrix}
-1&1&0\\
-1&0&1
\end{bmatrix}
$$










