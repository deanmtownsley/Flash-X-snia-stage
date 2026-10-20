!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant-diff/Conductivity
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
!!               real, intent(OUT)  :: isochoricCond,
!!               real, intent(OUT)  :: diff_coeff,
!!               integer,intent(IN) :: component)
!!
!! DESCRIPTION
!!   
!!   Applies conductivity with constant coefficient of diffusivity
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond  :   isochoric conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are reqested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!
!!
!!***

subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)
  
  use Conductivity_data,ONLY : cond_diffConstant
  use cond_commonData,ONLY   : cond_useConductivity
  use Eos_interface, ONLY : Eos
  use cond_interface,ONLY: cond_getCv
  implicit none
  
#include "constants.h"  
#include "Simulation.h"
#include "Eos.h"
  
  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, isochoricCond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  if (cond_useConductivity) then
     diff_coeff = cond_diffConstant
     isochoricCond = (cond_diffConstant * &
          cond_getCv(xden,xtemp,component,massFrac) * xden)

  else
     isochoricCond = 0.0
     diff_coeff = 0.0
  end if

end subroutine Conductivity
