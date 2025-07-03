!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/cond_commonInit
!!
!! NAME
!!
!!  cond_commonInit
!!
!! SYNOPSIS
!!
!!  call cond_commonInit()
!!
!! DESCRIPTION
!!
!!  Initialize shared data that can be used by all Conductivity implementations
!!
!! ARGUMENTS
!!
!!  none
!!
!!
!!***
subroutine cond_commonInit
  use cond_commonData, ONLY: cond_useConductivity, cond_meshMe, &
       cond_mele, cond_boltz, cond_hbar, cond_qele, &
       cond_navo, cond_c
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use Driver_interface, ONLY : Driver_getMype

  implicit none

#include "constants.h"

  call RuntimeParameters_get("useConductivity", cond_useConductivity)

  ! Everybody should know this
  call Driver_getMype(MESH_COMM,cond_meshMe)

  ! Set physical constants:
  call PhysicalConstants_get("electron mass",cond_mele)
  call PhysicalConstants_get("Boltzmann",cond_boltz)
  call PhysicalConstants_get("electron charge",cond_qele)
  call PhysicalConstants_get("Planck",cond_hbar)
  cond_hbar = cond_hbar/(2.0*PI)

  call PhysicalConstants_get("speed of light",cond_c)
  call PhysicalConstants_get("Avogadro", cond_navo)

  ! Any additional initialization in case anisotropic conductivity is used
  call cond_anisoInit()

end subroutine cond_commonInit
