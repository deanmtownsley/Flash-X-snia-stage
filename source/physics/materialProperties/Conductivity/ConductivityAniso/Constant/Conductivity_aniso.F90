!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/Constant/Conductivity_aniso
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
!!   This routine returns anisotropic conductivities
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
!!   If possible, the alternative routine Conductivity_anisoFullState should be called
!!   instead if this routine.
!!
!!   In this implementation, the value of component has no effect on the values
!!   returned in isochoricCond.
!!
!!  SEE ALSO
!!
!!   Conductivity_anisoFullState
!!***

subroutine Conductivity_aniso(xtemp,xden,massfrac,isochoricCond,diff_coeff,component)

  ! a generic conductivity routine.  Just return constant coefficients for anisotropic
  ! isochoric conductivity.
  ! Note that the elements of isochoricCond are constant but the elements ofdiff_coeff,
  ! containing coefficients of diffusivity, are not (in general)!
  
!!$  use Conductivity_data, ONLY: cond_constantIsochoric - the corresponding parameter for isotropic Constant conductivity
  use cond_anisoData,  ONLY:   cond_constantParallel, cond_constantPerpendicular, cond_constantCross
  use cond_commonData, ONLY:   cond_useConductivity
  use cond_interface,  ONLY: cond_getCv
  implicit none

#include "constants.h"
#include "Flash.h"

  real, intent(IN) :: xtemp, xden
  real, intent(OUT) ::  diff_coeff(3), isochoricCond(3)
  real, dimension(NSPECIES), intent(IN) :: massfrac
  integer,intent(IN) :: component

  real    :: cv
  
  if (cond_useConductivity) then

     isochoricCond(1) = cond_constantParallel
     isochoricCond(2) = cond_constantPerpendicular
     isochoricCond(3) = cond_constantCross

     cv            = cond_getCv(xden,xtemp,component,massFrac)
     diff_coeff(:) = isochoricCond(:) /(xden*cv)
  else
     isochoricCond(:) = 0.0
     diff_coeff(:) = 0.0
  end if
 

end subroutine Conductivity_aniso
