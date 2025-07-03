!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/Constant/Conductivity_anisoFullState
!!
!! NAME
!!  Conductivity_anisoFullState
!!
!! SYNOPSIS
!!  call Conductivity_anisoFullState(real(in)    :: solnVec(NUNK_VARS),
!!                     OPTIONAL,real(out)   :: isochoricCond(3),
!!                     OPTIONAL,real(out)   :: diffCoeff(3),
!!                     OPTIONAL,integer(in) :: component)
!!
!! DESCRIPTION
!!
!!  Returns thermal conductivity and/or diffusivity coefficients.
!!
!!  The stub implementation for Conductivity_anisoFullState returns isochoricCond = 0
!!  and diffCoeff = 0.
!!
!!  This implementation just calls the alternative interface Conductivity_aniso
!!  and returns its output if diffCoeff is requested.
!!
!!  The following applies in the context of the"Constant" implementation of
!!  anisotropic Conductivity:
!!  This implementation returns constant values for the three conductivity
!!  coefficients, which come directly from runtime parameters; these constant
!!  values are currently independent of the 'component' requested.
!!  Note that in this context diffusion coefficients are not necessarily constant,
!!  since they are related to conductivities by factors that involve the state
!!  of the material, in particular densities and specific heat coefficients,
!!  potentially in a component-specific way (particularly for 3T applications).
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   isochoricCond  :   isochoric conductivities (parallel, perpendicular, cross)
!!   diffCoeff :   diffusion coefficients ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!
!!
!!***

#include "Flash.h"

subroutine Conductivity_anisoFullState(solnVec,isochoricCond,diffCoeff,component)
  use cond_interface,  ONLY: cond_getCv
  use cond_anisoData,  ONLY:   cond_constantParallel, cond_constantPerpendicular, cond_constantCross
  use cond_commonData, ONLY:   cond_useConductivity

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff(3)
  real,    OPTIONAL, intent(OUT)  :: isochoricCond(3)
  integer, OPTIONAL, intent(IN) :: component

  real :: isochoricCondLoc(3), diffCoeffLoc(3)
  real :: cv

  isochoricCondLoc(:) = 0.0
  diffCoeffLoc(:) = 0.0
  
#ifdef DENS_VAR

  if (cond_useConductivity) then
     isochoricCondLoc(1) = cond_constantParallel
     isochoricCondLoc(2) = cond_constantPerpendicular
     isochoricCondLoc(3) = cond_constantCross
     if (present(diffCoeff)) then
        cv = cond_getCv(component=component,solnVec=solnVec)
        diffCoeffLoc(:) = isochoricCondLoc(:)/(solnVec(DENS_VAR)*cv)
     end if
  end if
#endif

  if(present(isochoricCond)) isochoricCond(:) = isochoricCondLoc(:)
  if(present(diffCoeff)) diffCoeff(:) = diffCoeffLoc(:)

end subroutine Conductivity_anisoFullState
