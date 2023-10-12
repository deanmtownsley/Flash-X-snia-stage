# Helmholtz + WeakLib Hybrid Equation of State

This hybrid implementation utilizes the WeakLib EOS at densities $\rho > \rho_\text{HI}$ and the Helmholtz EOS at densities $\rho < \rho_\text{LO}$.  In the transition region $\rho_\text{LO} \le \rho \le \rho_\text{HI}$, for a dependent variables $U_\text{WL}$ and $U_\text{HLM}$ obtained from Helmholtz and WeakLib, respectively, the hybrid EOS calculates the average:

$$
   U_\text{HYB} = \frac{c_\text{LO}U_\text{HLM} + c_\text{HI}U_\text{WL}}{c_\text{LO} + c_\text{HI}}
$$

where the weights are

$$
   c_\text{LO} = \left(\frac{\rho_\text{LO}}{\rho}\right)^n,\quad c_\text{HI} = \left(\frac{\rho}{\rho_\text{HI}}\right)^n
$$

s.t. $c_\text{LO},c_\text{HI} \le 1$, and $n > 1$.  Increasing the value of $n$ will increasingly favor the "closer" density boundary's EOS in the weighted average.

## Energy Offsets

The Helmholtz EOS requires the thermal