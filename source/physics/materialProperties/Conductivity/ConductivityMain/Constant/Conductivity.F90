!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant/Conductivity
!!
!! NAME
!!
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
!!   This routine implements constant conductivity.
!!
!!   This implementation returns a constant value for the conductivity
!!   coefficient, which comes directly from the runtime parameter
!!   cond_constantIsochoric. This value returned is currently
!!   independent of the 'component' requested.
!!
!!   Note that in this context the diffusion coefficient is not necessarily constant,
!!   since it is related to the conductivity by a factor that involves the state
!!   of the material, in particular densities and specific heat coefficients,
!!   potentially in a component-specific way (particularly for 3T applications).
!!
!! ARGUMENTS
!!
!!   xtemp      :   temperature (in K)
!!   xden       :   density (in g/cm**3)
!!   massfrac   :   mass fractions of the composition
!!   isochoricCond       :   isochoric conductivity
!!   diff_coeff :   diffusion coefficient ( = isochoricCond/(rho*cv))
!!   component  :   In 3T applications, select component for which conductivity
!!                  and diffusivity are requested, 1 for ions, 2 for electrons, 
!!                  3 for radiation.
!!
!!  NOTES
!!
!!   The returned isochoricCond is constant but the coefficient of
!!   diffusivity, diff_coeff, is not (in general).
!!***

subroutine Conductivity(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)

  ! a generic conductivity routine.  Just return a constant isochoric conductivity.
  ! But while isochoricCond is constant, the coefficient of diffusivity, diff_coeff,
  ! is not (in general)!
  
  use Conductivity_data, ONLY: cond_constantIsochoric, cond_useConductivity
  use cond_interface,    ONLY: cond_getCv
  implicit none

#include "Flash.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff, isochoricCond
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component
  
  if (cond_useConductivity) then

     isochoricCond = cond_constantIsochoric
     diff_coeff = isochoricCond/(xden*cond_getCv(xden,xtemp,component,massFrac))

  else
     isochoricCond = 0.0
     diff_coeff = 0.0
  end if
 

end subroutine Conductivity
