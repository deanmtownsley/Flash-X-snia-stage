!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/Constant/cond_anisoData
!!
!! NAME
!!
!!  cond_anisoData
!!
!! SYNOPSIS
!!
!!  use cond_anisoData
!!
!! DESCRIPTION
!!
!!  Stores additional local data for the Anisotropic Conductivity implementation.
!!
!!
!!***

Module cond_anisoData

  implicit none

  real, save :: cond_constantParallel
  real, save :: cond_constantPerpendicular
  real, save :: cond_constantCross
  
end Module cond_anisoData
