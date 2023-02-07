# Z4c Spacetime Solver

This unit implements a Z4c solver for Einstein's Equations of General Relativity

# TODO List

- [ ] Description of the solver in `README.md`
- [X] Add additional runtime parameters for gauge+slicing choices
- [X] Modify the RHS calculations in `Spacetime_molFastRHS_tile` to reflect the gauge+slicing choices (e.g. enforcing a zero shift)
- [X] Add conversion of Z4c variables to ADM variables
- [ ] Add matter (stress-energy tensor) source terms
- [ ] Add procedure for handling boundary conditions for Z4c vector and tensor variables that will be called from a `Simulation` unit's `Grid_bcApplyToRegionSpecialized`
- [ ] Add procedure for handling Z4c-specific refinement criteria that will be called from a `Simulation` unit's `Grid_markRefineSpecialized`
- Update code-generator for more readable output?
