!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw/Conductivity
!!
!! NAME
!!
!!  Conductivity
!!
!! SYNOPSIS
!!
!!  Conductivity(real, intent(IN)  :: xtemp,
!!               real, intent(IN)  :: xden,
!!               real, dimension(NSPECIES), intent(IN)  :: massfrac,
!!               real, intent(OUT)  :: cond,
!!               real, intent(OUT)  :: diff_coeff,
!!               integer,intent(IN)     :: component)
!!
!! DESCRIPTION
!!
!!   Thermal conductivity and diffusion coefficient given by a power law;
!!   this includes Spitzer (1962)
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
!!  NOTES
!!   See: Spitzer L. 1962, In: `Physics of fully ionized gases', 
!!   (New York: Wiley Interscience)
!!
!!  MODIFICATION HISTORY
!!
!!   Initial: K. Weide March 2010
!!            K. Weide Feb   2019 - simplified, usecond_getCv
!!
!!***



subroutine Conductivity(xtemp,xden,massfrac,cond,diff_coeff,component)

  use Conductivity_data, ONLY: cond_useConductivity,  &
       cond_TemperatureExponent, cond_K0, cond_alpha, &
       cond_DensityExponent 
  use cond_interface,ONLY: cond_getCv
  implicit none
  
#include "Flash.h"
  
  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, cond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real :: ck
  real :: dexp
  real :: cexp

  if (cond_useConductivity) then
     cexp = cond_TemperatureExponent
     dexp = cond_DensityExponent
     ck   = cond_K0
     if (xtemp*xden .NE. 0.0) then
        cond       = ck*(xtemp**cexp)*(xden**dexp)
        diff_coeff = cond/(xden*cond_getCv(xden,xtemp,component,massFrac))
     else
        cond = 0.0
        diff_coeff = 0.0
     end if
  
  else
     cond = 0.0
     diff_coeff = 0.0
  end if

end subroutine Conductivity
