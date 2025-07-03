!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/cond_commonData
!!
!! NAME
!!
!!  cond_commonData
!!
!! SYNOPSIS
!!
!!  use cond_commonData
!!
!! DESCRIPTION
!!
!!  Private local data that is common to all Conductivity implementations
!!
!!
!!***

Module cond_commonData

  implicit none

  logical, save::  cond_useConductivity
  
  integer, save::  cond_meshMe

  ! Storage for some physical constants:
  real, save :: cond_mele  ! Electron mass (grams)
  real, save :: cond_boltz ! Boltzmann constant (ergs/K)
  real, save :: cond_hbar  ! Planck's constant over 2 PI (erg*s)
  real, save :: cond_qele  ! Proton charge (esu)

  real, save :: cond_navo  ! Avogadros number
  real, save :: cond_c     ! speed of light (cm/s)
  
end Module cond_commonData
