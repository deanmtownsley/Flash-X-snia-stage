# Brusselator Test

## Setup line

```bash
./setup Brusselator -auto +ug +nofbs -1d +MoLMR
```

The Brusselator test (see [Chinomona & Reynolds, 2021](https://arxiv.org/abs/2007.09776)) is a stiff advection-diffusion-reaction problem intended to demonstrate the capabilities of the method-of-lines (MoL) multi-rate time integrator.  This simulation evolves the system of PDEs

$$
\begin{align*}
\frac{\partial u}{\partial t} &= \rho\frac{\partial u}{\partial x} + \alpha\frac{\partial^2 u}{\partial x^2} + \left[a - (w-1)u + u^2v\right]\\\\
\frac{\partial v}{\partial t} &= \rho\frac{\partial v}{\partial x} + \alpha\frac{\partial^2 v}{\partial x^2} + \left[wu - u^2v\right]\\\\
\frac{\partial w}{\partial t} &= \rho\frac{\partial w}{\partial x} + \alpha\frac{\partial^2 w}{\partial x^2} + \left[\frac{b-w}{\epsilon}-wu\right]
\end{align*}
$$

on the domain $x\in[0,1]$, with stationary boundary conditions

$$
\frac{\partial u}{\partial t}\biggr\rvert_{x=0,1} = \frac{\partial v}{\partial t}\biggr\rvert_{x=0,1} = \frac{\partial w}{\partial t}\biggr\rvert_{x=0,1} = 0
$$

and initial conditions

$$
\begin{align*}
u(0,x) &= a + 0.1\sin(\pi x)\\
v(0,x) &= b/a + 0.1\sin(\pi x)\\
w(0,x) &= b + 0.1\sin(\pi x)
\end{align*}
$$

- The parameters $a$, $b$, $\alpha$, $\rho$, $\epsilon$ are all runtime inputs (see the included `flash.par`)
- The advective terms (first term) are treated explicitly during the "slow" integration steps, and computed via upwind-differencing
- The diffusive terms (second term) are treated implicitly during the "slow" integration steps, and computed with 2nd-order centered-differennces with the resulting tri-diagonal matrix equation solved on the grid with the `dgtsv` routine in LAPACK
- The reaction terms (third term in brackets) are treated explicitly over a serious of "fast" integration steps
- The time step is controlled as $\Delta t = 0.1\cdot 2^{-k}$, where $k$ is configurable at runtime
- The first and last cells are configured to be centered at $x = 0,1$, and the stationary boundary conditions are enforced during calculation of the RHS terms.
