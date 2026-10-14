!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/SpitzerHighZ/Conductivity_data
!!
!! NAME
!!
!!  Conductivity_data
!!
!! SYNOPSIS
!!
!!  use Conductivity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for the Spitzer Conductivity implementation.
!!
!!
!!***

Module Conductivity_data

  use cond_commonData, ONLY: cond_useConductivity, cond_meshMe
  use cond_commonData, ONLY: cond_mele, cond_boltz, cond_hbar, cond_navo, cond_qele
  !includes the following:
  ! real, save :: cond_mele  ! Electron mass (grams)
  ! real, save :: cond_boltz ! Boltzmann constant (ergs/K)
  ! real, save :: cond_hbar  ! Planck's constant over 2 PI (erg*s)
  ! real, save :: cond_qele  ! Proton charge (esu)

  implicit none

  ! Storage for additional constants:
  
end Module Conductivity_data
