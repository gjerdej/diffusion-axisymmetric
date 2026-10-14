# Cylindrical Axisymmetric Diffusion in PRISMS

We wish to solve the diffusion equation

$$\frac{\partial c}{\partial t} = - \nabla \cdot (D\nabla c)$$

With backward Euler time stepping

$$c^{n + 1} - c^{n} = - \Delta t\nabla \cdot \left( D\nabla c^{n + 1} \right)$$

If we multiply by the test function $w$ and integrate over the cylindrical domain

$$\int_{}^{}{w\left( c^{n + 1} - c^{n} \right)rdrd\theta dz} = - \int_{}^{}{\Delta tw\nabla \cdot \left( D\nabla c^{n + 1} \right)rdrdz}$$

If the problem is axisymmetric, we can perform the $\theta$ integration. The result is a common factor of $2\pi$ that may be cancelled. Integrating by parts and assuming the boundary term is zero

$$\int_{}^{}{(wc^{n + 1} - \nabla w \cdot \left( D\nabla c^{n + 1} \right))rdrdz} = \int_{}^{}{wc^{n}rdrdz}$$

## Initial and Boundary Conditions for Verification

We set all boundary conditions to natural and apply a non-uniform initial condition. Eventually, the solution will be uniform everywhere, since no flux is entering or leaving. Choosing the initial condition carefully will allow us to find an analytical solution to compare with our numerical solution.

The boundary conditions permit a separable solution (i.e. $c = \rho(r)\zeta(z)\tau(t)$ is a valid solution). The general form of the solution is

$$c(r,z,t) = \sum_{n = 0,m = 0}^{\infty}{a_{mn}J_{0}\left( k_{n}r \right)\left( B_{m}\cos\left( l_{m}z \right) + C_{m}\sin\left( l_{m}z \right) \right)\exp\left( - D\left( k_{n}^{2} + l_{m}^{2} \right)t \right)}$$

where $J_{0}(x)$ is the Bessel function of the first kind and zeroth order. If we apply the remaining boundary conditions

$$c(r,z,t) = \sum_{n = 0,m = 0}^{\infty}{a_{mn}J_{0}\left( \frac{\delta_{n}r}{R} \right)\sin\left( \frac{m\pi z}{H} \right)\exp\left( - D\left( \left( \frac{\delta_{n}}{R} \right)^{2} + \left( \frac{m\pi}{H} \right)^{2} \right)t \right)}$$

where $\delta_{n}$ are the zeros of $J_{0}(x)$, with $\delta_{0} = 0$.

We can choose our initial condition such that only the constant mode and one non-constant mode are nonzero

$$c(r,z,0) = c_{0} + AJ_{0}\left( \frac{\delta_{1}r}{R} \right)\sin\left( \frac{\pi z}{H} \right)$$

By orthogonality, $a_{00} = c_{0}$, $a_{11} = A$, and all other $a_{nm}$ are zero.

$$c(r,z,0) = c_{0} + AJ_{0}\left( \frac{\delta_{0}r}{R} \right)\sin\left( \frac{\pi z}{H} \right)\exp\left( - D\left( \left( \frac{\delta_{n}}{R} \right)^{2} + \left( \frac{m\pi}{H} \right)^{2} \right)t \right)$$

The analytical solution is included as a post-processed field.
