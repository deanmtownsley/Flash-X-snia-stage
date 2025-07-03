!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity
!!
!! NAME
!!
!!  Conductivity
!!
!! SYNOPSIS
!!
!!  call Conductivity(real, intent(IN)  :: xtemp,
!!                    real, intent(IN)  :: xden,
!!                    real, dimension(NSPECIES), intent(IN)  :: massfrac,
!!                    real, intent(OUT)  :: cond,
!!                    real, intent(OUT)  :: diff_coeff,
!!                    integer(in) :: component)
!!
!! DESCRIPTION
!!
!! Computes the Spitzer electron conductivity for all materials,
!! including those with Z > 1. The specific equations used here all
!! come from "The Physics of Inertial Fusion" by Atzeni.
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   cond       :   conductivity
!!   diff_coeff :   diffusion coefficient ( = cond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are reqested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!***

#include "Flash.h"

subroutine Conductivity(xtemp,xden,massfrac,cond,diff_coeff, component)
  use Conductivity_data, ONLY: cond_useConductivity, &
       cond_mele, cond_boltz, cond_qele, cond_navo
  use Eos_interface, ONLY: Eos_getAbarZbar
  use cond_interface,ONLY: cond_getCv
  implicit none
  
#include "constants.h"  
#include "Eos.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, cond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real :: nele, ll
  real :: abar, zbar

  real, parameter :: cexp = 2.5

  if (cond_useConductivity) then
     call Eos_getAbarZbar(abar=abar,zbar=zbar,massFrac=massfrac)

     nele = zbar * xden * cond_navo / abar
     call cond_logLambda(xtemp, nele, xtemp, zbar, ll)
     
     cond = (8.0/PI)**1.5*cond_boltz**3.5 / (sqrt(cond_mele)*cond_qele**4) * &
          xtemp**cexp / (ll * (zbar + 3.3))

     diff_coeff = cond/(xden*cond_getCv(xden,xtemp,component,massFrac))
  
  else
     cond = 0.0
     diff_coeff = 0.0
  end if

end subroutine Conductivity
