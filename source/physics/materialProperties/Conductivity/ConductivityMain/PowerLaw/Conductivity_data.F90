!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw/Conductivity_data
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
!!  Stores the local data for the power law conductivity implementation.
!!
!!
!!***
Module Conductivity_data

  use cond_commonData, ONLY: cond_useConductivity, cond_meshMe

  implicit none

  real, save :: cond_TemperatureExponent, cond_K0, cond_alpha
  real, save :: cond_DensityExponent

end Module Conductivity_data
