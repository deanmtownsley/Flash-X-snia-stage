!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw-gray/Conductivity_init
!!
!! NAME
!!
!!  Conductivity_init
!!
!! SYNOPSIS
!!
!!  Conductivity_init()
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!
!!***

#include "constants.h"

subroutine Conductivity_init()

  use Conductivity_data, ONLY: cond_TemperatureExponent, cond_K0, cond_alpha, Raddiff_K0r, Raddiff_TemperatureExponent
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  implicit none

  

  call cond_commonInit()

  call RuntimeParameters_get("cond_TemperatureExponent", cond_TemperatureExponent)
  call RuntimeParameters_get("cond_K0", cond_K0)
 
  call RuntimeParameters_get("Raddiff_K0r", Raddiff_K0r)
  call RuntimeParameters_get("Raddiff_TemperatureExponent", Raddiff_TemperatureExponent)

  

  cond_alpha = cond_K0


end subroutine Conductivity_init

