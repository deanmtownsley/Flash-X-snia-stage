!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/Conductivity_aniso
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
!! Computes the Spitzer electron conductivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!! The true anisotropic conductivities from Braginskii are calculated
!! in Conductivity_fullState.F90 
!! 
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond: conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!  NOTES
!!
!!   Currently, the alternative routine Conductivity_anisoFullState has to be called
!!   instead if this routine to actually get correct anisotropic coefficients.
!!
!!   In this implementation, the value of component has no effect on the values
!!   returned in isochoricCond.
!!
!!  SEE ALSO
!!
!!   Conductivity_anisoFullState
!!***

#include "Simulation.h"

subroutine Conductivity_aniso(xtemp,xden,massfrac,isochoricCond,diff_coeff, component)
  use cond_commonData, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele
  use cond_commonData, ONLY: cond_navo
  use Eos_interface, ONLY: Eos_getAbarZbar
  use cond_interface,ONLY: cond_getCv
  implicit none
  
#include "constants.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff(3), isochoricCond(3)
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real :: nele, ll
  real :: abar, zbar
  real :: cv

  real, parameter :: cexp = 2.5

  if (cond_useConductivity) then
     call Eos_getAbarZbar(abar=abar,zbar=zbar,massFrac=massfrac)

     nele = zbar * xden * cond_navo / abar
     call cond_logLambda(xtemp, nele, xtemp, zbar, ll)
     
     isochoricCond(:) = (8.0/PI)**1.5*cond_boltz**3.5 / (sqrt(cond_mele)*cond_qele**4) * &
          xtemp**cexp / (ll * (zbar + 3.3))

     cv = cond_getCv(xden,xtemp,component,massFrac)
     diff_coeff(:) = isochoricCond(:)/(xden*cv)
  
  else
     isochoricCond(:) = 0.0
     diff_coeff(:) = 0.0
  end if

end subroutine Conductivity_aniso
