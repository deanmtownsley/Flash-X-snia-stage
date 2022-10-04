# Method of Lines Time Integration Unit

The purpose of this README is to give a very general overview of the current state of `MoL` and record some of the outstanding issues/comments/plans related to discussions in the PR and during various related meetings.

This unit provides a method of lines (`MoL`) time integrator as an alternative to the operator
split scheme available by-default in `Driver`.  The integrators contained in `MoL` allow for tightly
coupled equations split up amongst one-or-more physics units to be evolved consistently at 
a per-stage level during a time step.

Currently, `MoL` provides three distinct integrators:
- `ERK`: explicit Runge-Kutta methods
- `FBE`: forward-backward Euler implicit-explicit (IMEX) method
   - Note: this may be extended to a more general IMEX integrator to include higher-order methods
- `MR`: multi-rate methods that couple an IMEX method on a "slower" time scale with an explicit method (currently - possibly implicit or IMEX in the future if needed) on a "faster" time scale

## Usage Notes

Using one of the specific types of integrators above requires either:
- A `REQUIRES numericalTools/MoL/MoLMain/[ERK|FBE|MR]` in a simulation's (or physics unit's) `Config` file
- A `--with-unit=numericalTools/MoL/MoLMain/[ERK|FBE|MR]` in the setup line

When one of these are present, a `MoL`-specific set of `Driver` will also be loaded, including a modified `Driver_evolveAll` and a set of procedures `dr_mol*` that aggregate various calls to physics units to perform:
- Calculation of right-hand side (RHS) terms of the equations being evolved
- Post-update (per-stage and per-timestep) work, such as `con2prim` operation
- Specialized updates that are not included more generally in a specific integrator, e.g. an implicit update of a stiff source term (these are currently kicked back to the appropriate physics unit(s) to to exploit any symmetries that may exist in the equations)

Currently, there are three sets of interfaces in `Simulation`, `Hydro`, and `RadTrans` (of the form `Unit_mol*`) that provide a simple means of connecting a physics unit to `MoL`.

Additionally, usage of `MoL` will require the use of `--with-unofficial=numericalTools/MoL` on the setup line until the `UNOFFICIAL` tag is removed from the `MoL` unit.

## Evolved Variables

Evolved variables can be any declared variable, i.e. `VARIABLE <NAME> ...` in a `Config` file.  For now, there is a start-up requirement of informing `MoL` of which variables are evolved and obtaining a corresponding index for each of these variables in the `MoL`-specific memory data structures.  For example, a variable for the x-component of a momentum density can be "registered" with `MoL` by
```fortran
call MoL_registerVariable("momx", MOMX_VAR, MOMX_RHS)
```
where `MOMX_VAR` is a setup-time generated index for the variable in `Simulation.h`, and `MOMX_RHS` is an integer variable that will be set to the index of the variable in MoL's RHS memory structures.

### _Future Changes to Evolved Variables_

The need for variable registration is currently necessary due to how `MoL` manages its own scratch-memory for various intermediate states that are required in some integration methods.  This is not an ideal solution, as the presence and knowledge of the evolved variables are known at setup-time.  It would be best to remove this level of redundancy between setup-time configurations and runtime startup.  One potential change may include:

- Addition of a configuration syntax of the form `EVOLVEDVAR <NAME>` for declaring evolved variables in a `Config` file
- For every declared evolved variable, the setup process will generate a unique identifer in `UNK` of a form `NAME_VAR` (or `NAME_EVOL` perhaps) or possibly a physical data group specifically for evolved variables
- For every declared evolved variable, the setup process will generate a series of unique identifers in `Grid`-managed scratch memory (setup-time or allocatable) of a form `NAME_RHS_X` where `X` can be a number/identifier for the integration stage and/or type of RHS that this corresponds to (explicit, implicit, etc.)
- The setup process will be able to count the total number of evolved variables and other necessary information that an integrator may require
- For ease-of-use and reducing redundancies in the code about the details of any particular integration method, a set of setup-time defines (possibly for each variable) can be made to inform the various physics units how many RHS variables exist for all (or each) variable, allowing for a compile-time header approach in each unit to generate a indexable-list of the RHS variable identifers, e.g. using `MOMX` as an example this may look like
```c
#define MOMX_RHS_BEGIN MOMX_RHS_1
#define MOMX_RHS_END   (MOMX_RHS_BEGIN + MOL_NRHS - 1)
```
and
```fortran
! in some data module
integer, dimension(MOL_NRHS), save :: MOMX_RHS_INDS
...
! in some init procedure
MOMX_RHS_INDS = [(i, i = MOMX_RHS_BEGIN, MOMX_RHS_END)]
```
- An "active RHS" identifier in the range `[1, MOL_NRHS]` could be passed to the procedures responsible for calculating the RHS terms and used like:
```fortran
subroutine Hydro_molExplicitRHS(t, activeRHS)
   use MoL_interface, only: MoL_getDataPtr, ...
   use hy_data, only: MOMX_RHS_INDS, ...
   ...

   implicit none

   real, intent(in) :: t
   integer, intent(in) :: activeRHS

   integer :: MOMX_RHS
   real, pointer :: rhs(:,:,:)
   ...

   MOMX_RHS = MOMX_RHS_INDS(activeRHS)
   ...

   ! Somewhere in a tile-loop
   call MoL_getDataPtr(tileDesc, rhs, activeRHS)
   ...

   ! Somewhere in an i/j/k loop
   rhs(MOMX_RHS, i, j, k) = rhs(MOMX_RHS, i, j, k) + ...

   ...
end subroutine Hydro_molExplicitRHS
```
- There are times when some level of abtraction as to the type and stage that `activeRHS` corresponds to is beneficial, e.g. if a physics unit implements a "fast" RHS procedure, but if the ERK integrator is selected, there will be no specific "fast" RHS storage allocated, so this procedure will be called on to add its "fast" contribution to the same location as the explicit RHS terms
- This type of evolved and RHS variable indexing could also simplify some of the internal `MoL` operations on these variables (e.g. linear combinations of intermediate state RHS terms) by removing the explicit need for indexing mapping between `UNK` and the current version of MoL's scratch memory

### Scratch Memory in MoL

`MoL` currently implements a variation of `GridAllocatableScratches` using a different rank-layout that is better suited for `MoL`'s iteration and combinations of its intermediate-state memory structures.  This is meant only as a temporary solution, as there is little-to-no benefit having to replicate a subset of `Grid`'s functionality for each and every different backend.  A few possible ideas on how to best replace this include:

- Extend `GridAllocatableScratches` to included an additional layout for use in `MoL`
- Utilize setup-time `SCRATCHVAR`s via the proposed (or similar to) configuration syntax above
- Potentially utilize multiple physical data groups alongside some (smaller) additions to the configuration syntax to generate and/or group evolved and RHS variables

Some of the clunkiness of the current implementation (e.g. the use of the older style `getDataPtr` call that is not part of the tile descriptor) is intentional as it will allow for fewer/smaller changes to be made to code outside of `MoLMemory` once a more-suitable solution for the scratch memory requirements is brought into `MoL`.

## Other Notes and Action Items
I am relatively certain I have missed adding in some open topics from discussion in/about the PR, and I will continue to add/update these as necessary.

I plan on generating a more user-oriented guide on using `MoL` in a physics unit as we begin working on updating and/or adding units that utilize `MoL` (the in-development spacetime solver for example).  Currently, the simulation `Simulation/SimulationMain/Brusselator` demonstrates how to utilize `MoL` for a simple advection-diffusion-reaction problem.