!!****f* source/physics/materialProperties/Conductivity/Conductivity_anisoFullState
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
  use Conductivity_interface, ONLY: Conductivity_aniso

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff(3)
  real,    OPTIONAL, intent(OUT)  :: isochoricCond(3)
  integer, OPTIONAL, intent(IN) :: component

  real :: isochoricCondLoc(3), diffCoeffLoc(3)

  isochoricCondLoc(:) = 0.0
  diffCoeffLoc(:) = 0.0

  if(present(isochoricCond)) isochoricCond(:) = isochoricCondLoc(:)
  if(present(diffCoeff)) diffCoeff(:) = diffCoeffLoc(:)

end subroutine Conductivity_anisoFullState
