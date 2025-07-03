!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/PowerLaw-gray/Conductivity_data
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
!!  Stores the local data for the Powerlaw gray Conductivity implementation.
!!
!!
!!***
Module Conductivity_data

  use cond_commonData, ONLY: cond_useConductivity, cond_meshMe
  !see diff_* in Diffuse_data for  useRaddiffusion, useElecond
  

  implicit none

  real, save :: cond_TemperatureExponent, cond_K0, cond_alpha, Raddiff_K0r, Raddiff_TemperatureExponent

end Module Conductivity_data
