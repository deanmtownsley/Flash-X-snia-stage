# Helmholtz + WeakLib Hybrid Equation of State

This hybrid implementation utilizes the WeakLib EOS at densities $\rho > \rho_\text{HI}$ and the Helmholtz EOS at densities $\rho < \rho_\text{LO}$.  In the transition region $\rho_\text{LO} \le \rho \le \rho_\text{HI}$, for a dependent variables $U_\text{WL}$ and $U_\text{HLM}$ obtained from Helmholtz and WeakLib, respectively, the hybrid EOS calculates the average:

$$
   U_\text{HYB} = (1-c)U_\text{HLM} + c\,U_\text{WL}
$$

where the weight $c$ is

$$
   c = \frac{\rho - \rho_\text{LO}}{\rho_\text{HI} - \rho_\text{LO}}
$$

## Energy Offsets

The Helmholtz EOS requires the thermal energy $E_\text{th}$ as input, so the internal energy $E_\text{int}$ must be offset by the rest mass energy and any other shifts associated with the nuclear EOS tables.  For WeakLib, there is a shift of $\Delta E_\text{WL} = 8.9$ MeV/nucleon.  The energy shift due to the rest mass energy is

$$
   \Delta E_\text{M} = -\left[Y_e(m_n - m_p - m_e)c^2 + \sum_i Y_i B_i\right]
$$

where $Y_i$ and $B_i$ are the abundance and binding energy of the $i$-th evolved/tracked species.  The thermal energy is then

$$
   E_\text{th} = E_\text{int} - \Delta E_\text{M} - \Delta E_\text{WL}
$$
