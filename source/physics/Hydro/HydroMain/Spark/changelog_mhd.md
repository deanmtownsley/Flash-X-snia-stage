# Changes to tvd2 hy_recon_riemann.ini

- Specified ```Bn_glm``` and ```Psi_glm``` for GLM 2x2 sub-system procedure
- Added the following changes:

```Fortran
Bn_glm  = 0.5*(VL(HY_MAGX+dir-1)+VR(HY_MAGX+dir-1)) - 0.5/C_hyp*(VR(HY_PSIB) - VL(HY_PSIB))
Psi_glm = 0.5*(VL(HY_PSIB)+VR(HY_PSIB)) - 0.5*C_hyp*(VR(HY_MAGX+dir-1)-VL(HY_MAGX+dir-1))
VL(HY_PSIB) = Psi_glm
VR(HY_PSIB) = Psi_glm
velNL = VL(HY_VELX+dir-1)
velNR = VR(HY_VELX+dir-1)
```

```Fortran
#ifdef SPARK_GLM
 ! The exact fluxes for the 2x2 GLM sub-system
 Fstar(HY_FMGX+dir-1) = Psi_glm
 Fstar(HY_FPSI) = C_hyp*C_hyp*Bn_glm
#endif
```

# Changes to WENO hy_recon_riemann.ini

- Added the following changes:

```Fortran
VL(HY_PSIB) = Psi_glm
VR(HY_PSIB) = Psi_glm
velNL = VL(HY_VELX+dir-1)
velNR = VR(HY_VELX+dir-1)
```

# Experimental changes with WENO (NOT INCLUDED HERE)

- Note that abs_betadiff/(betaWeno(1:3) + epsilon)^n in hy_recon_riemann.ini is hard-coded with n=1, which is fine for MHD.
- However, we find that n=2 works better for pure hydro. This may be worth investigating later on.
