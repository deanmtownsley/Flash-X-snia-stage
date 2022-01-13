# Design Document for Spark Hydro solver

## Desired features:

- Runge-Kutta time integration of arbitrary order
- High-order spatial interpolation (with compact stencils)
- Fully conservative momentum equation with self-gravity
- Optional fully conservative energy equation with self-gravity or Mueller trick
- Species/mass scalar advection
- Pressure incorporated into fluxes
- Arbitrary geometry
- (dens*eps) reconstruction for face state thermodynamics
- CMA for species (not sure if it really helps...)
- We will stick with primitive variable reconstruction
- Ultimately, cell-center vector potential CT MHD
- Would like to address fine-coarse interpolation
- May need to consider "well balancing" but I think conservative gravity is the way to go... Likely reduces to same thing.
- Note that the average fluxes for multi-stage RK methods can be computed easily, but not exactly for low-storage methods...
- Not sure if we want to worry about the internal energy equation for eint << ekin condtions...
- For Riemann solvers, have top level routine hy_rk_riemann and #includes the actual solvers AND the prim2con prim2flx
- Would be nice to have genuine adaptive time stepping by computing different accuracy solutions.
- On gravity coupling: Can use the dphit/dt to estimate the GPOT for the intermediate stages
- time sub-cycling: can use the fact that solnData guardcell info does not change. Maybe?

## Some conventions:

- "+" and "-" are always relative to cell centers, "L" and "R" are always relative to cell interfaces
- files thatbegin with "hy_rk_" are major routines that are not contained by any other routine, have an interface in hy_rk_interface.F90, and are included in the Makefile. Files that are not prepended by "hy_rk_" are routines that are contained within one of the major routines.

## Variables to reconstruct:

- DENS
- VELX
- VELY
- VELZ
- PRES
- GAMC
- DENS*EINT
- GRAV(XYZ)
- SPECS
- MSCALARS

## Implementation plan:

1. Start with second-order TVD, no gravity, no mscalars [check]

  - HLLE/HLLC Riemann solver(s)
  - verify with standard tests (Sod, etc.)

2. Incorporate geometry support [check]

3. Extend to third-order LIMO3 interpolation [check]

  - re-run all tests

4. Extend to RK3 [check]

  - re-run tests

5. Incorporate mass scalar advection [check]

6. Incorporate gravity via source terms [check]

7. Shock detection [check]

8. Incorporate gravity with conservative momentum and energy equations to first-order

  - Begin gravity tests, compare to Kappeli and Mishra
  - This will not work for non-Poisson gravity...

9. High-order in time gravity (a la Jiang et al.)

  - recompute potential after each RK stage
  - maybe. this may not really be needed. check energy conservation

10. Adaptive time stepping using the energy equation evolution

11. Add CMA, compare to standard to see if any benefit

12. Gonna need some flux correction for production runs...

  - could treat the second-order gravity corrections along with flux corrections NO
  - Step 1: store summation of weighted fluxes from RK intermediate steps. This may require modification to the Grid routines saving the fluxes. Perhaps an optional flag?
  - Step 2: Modify Grid_conserveFluxes to return the DIFFERENCES between old and new fluxes.
  - Step 3: May need a new Grid routine to flush the fluxes from the Grid data structures
  - Step 4: Correct the block-edge zones with conservative fluxes. Easiest thing to do is make sure all the fluxes except the edge fluxes are zero then update the whole damn block... That's a lot of extra compute... Will just make a separate update routine...
