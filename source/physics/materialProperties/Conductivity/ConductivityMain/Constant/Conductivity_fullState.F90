!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant/Conductivity_fullState
!!
!! NAME
!!  Conductivity_fullState
!!
!! SYNOPSIS
!!  call Conductivity_fullState(real(in)    :: solnVec(NUNK_VARS),
!!                     OPTIONAL,real(out)   :: isochoricCond,
!!                     OPTIONAL,real(out)   :: diffCoeff,
!!                     OPTIONAL,integer(in) :: component)
!!
!! DESCRIPTION
!!
!!  Returns thermal conductivity and/or diffusivity coefficients.
!!
!!  The stub implementation for Conductivity_fullState returns isochoricCond = 0
!!  and diffCoeff = 0.
!!
!!  This implementation calls the alternative specific subroutine
!!  Conductivity() and returns its output if diffCoeff is present.
!!
!! ARGUMENTS
!!
!!   solnVec  :   solution state, a vector from UNK with all variables
!!   isochoricCond  :   isochoric conductivity
!!   diffCoeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!
!!
!!***

#include "Flash.h"

subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_data, ONLY: cond_constantIsochoric, cond_useConductivity
  use cond_interface,    ONLY: cond_getCv

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff
  real,    OPTIONAL, intent(OUT)  :: isochoricCond
  integer, OPTIONAL, intent(IN) :: component

  real :: isochoricCondLoc, diffCoeffLoc

  isochoricCondLoc = 0.0
  diffCoeffLoc = 0.0
  
#if defined(DENS_VAR) && defined(TEMP_VAR)
  if (cond_useConductivity) then
     isochoricCondLoc = cond_constantIsochoric
     if (present(diffCoeff)) then
        diffCoeffLoc = isochoricCondLoc/(solnVec(DENS_VAR) * &
             cond_getCv(component=component,solnVec=solnVec))
     end if
  end if
#endif

  if(present(isochoricCond)) isochoricCond = isochoricCondLoc
  if(present(diffCoeff)) diffCoeff = diffCoeffLoc

end subroutine Conductivity_fullState

