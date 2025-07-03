!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/Constant/Conductivity_data
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
!!  Stores the local data for the constant conductivity implementation.
!!
!!
!!***
Module Conductivity_data

  use cond_commonData, ONLY: cond_useConductivity, cond_meshMe

  implicit none
  
  real, save :: cond_constantIsochoric

end Module Conductivity_data
