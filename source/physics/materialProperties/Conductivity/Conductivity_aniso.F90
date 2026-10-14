!!****f* source/physics/materialProperties/Conductivity/Conductivity_aniso
!!
!! NAME
!!
!!  Conductivity_aniso
!!
!! SYNOPSIS
!!
!!  call Conductivity_aniso(real(IN)   :: xtemp,
!!                          real(IN)   :: xden,
!!                          real(IN)   :: massfrac(NSPECIES),
!!                          real(OUT)  :: isochoricCond(3),
!!                          real(OUT)  :: diff_coeff(3),
!!                          integer(in) :: component)
!!
!! DESCRIPTION
!!
!!  Returns thermal conductivity and diffusivity coefficients.
!!
!!  The stub implementation for Conductivity_aniso returns isochoricCond = 0
!!  and diff_coeff = 0.
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond  :   isochoric conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!
!!
!!***

subroutine Conductivity_aniso(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)
!! true stub

#include "Simulation.h"
  
  implicit none

  real, intent(IN)   :: xtemp
  real, intent(IN)   :: xden
  real, intent(IN)   :: massfrac(NSPECIES)
  real, intent(OUT)  :: diff_coeff(3)
  real, intent(OUT)  :: isochoricCond(3)
  integer,intent(IN) :: component


  !dummy values assigned in stub
  diff_coeff = 0.
  isochoricCond = 0.


  return 
end subroutine Conductivity_aniso

subroutine Conductivity_anisoFullState(solnVec,isochoricCond,diffCoeff,component)
  use Conductivity_interface, ONLY: Conductivity

  implicit none

  real,              intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff(3)
  real,    OPTIONAL, intent(OUT)  :: isochoricCond(3)
  integer, OPTIONAL, intent(IN) :: component

  if(present(isochoricCond)) isochoricCond(:) = 0.0
  if(present(diffCoeff)) diffCoeff(:) = 0.0

end subroutine Conductivity_anisoFullState

