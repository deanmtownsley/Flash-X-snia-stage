!!****f* source/physics/materialProperties/Conductivity/Conductivity
!!
!! NAME
!!  Conductivity
!!
!! SYNOPSIS
!!  call Conductivity(real(in)    :: xtemp,
!!                    real(in)    :: xden,
!!                    real(in)    :: massfrac(NSPECIES),
!!                    real(out)   :: isochoricCond,
!!                    real(out)   :: diff_coeff,
!!                    integer(in) :: component)
!!
!! DESCRIPTION
!!
!!  Returns thermal conductivity and diffusivity coefficients.
!!
!!  The stub implementation for Conductivity returns isochoricCond = 0
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

subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)
!! true stub

#include "Simulation.h"
  
  implicit none

  real, intent(IN)   :: xtemp
  real, intent(IN)   :: xden
  real, intent(IN)   :: massfrac(NSPECIES)
  real, intent(OUT)  :: diff_coeff
  real, intent(OUT)  :: isochoricCond
  integer,intent(IN) :: component


  !dummy values assigned in stub
  diff_coeff = 0.
  isochoricCond = 0.


  return 
end subroutine Conductivity

subroutine Conductivity_fullState(solnVec,isochoricCond,diffCoeff,component)
!! true stub

  implicit none

  real,    target,   intent(IN) :: solnVec(NUNK_VARS)
  real,    OPTIONAL, intent(OUT)  :: diffCoeff
  real,    OPTIONAL, intent(OUT)  :: isochoricCond
  integer, OPTIONAL, intent(IN) :: component

  if(present(isochoricCond)) isochoricCond = 0.0
  if(present(diffCoeff)) diffCoeff = 0.0

end subroutine Conductivity_fullState

