!!****ih* source/physics/Eos/EosMain/Hybrid/eos_hybridData
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!! NAME
!!
!!  eos_hybridData
!!
!!
!! SYNOPSIS
!!
!! use eos_hybridData
!!
!! DESCRIPTION
!!
!!  This is the data module for the hybrid (Helmholtz+WeakLib) Eos implementation.
!!  It stores all the runtime parameters, and all the unit scope
!!  data. Some of the unit scope data is fecthed by the wrapper layer
!!  from elsewhere in the code and some is local unit data common
!!  multiple functions in the unit
!!
!! PARAMETERS
!!
!!   These are the runtime parameters used by the hybrid implementation of the Eos unit.
!!
!!***

module eos_hybridData

   implicit none

   real, save :: eos_hybTransitionDensHi
   real, save :: eos_hybTransitionDensLo

   ! These were in old hybridData
   real, parameter :: amu_cgs = 1.66053873e-24
   real, parameter :: kb_erg = 1.380658e-16

   real, parameter :: ergPerKg_to_kbPerBaryon = amu_cgs/kb_erg
   real, parameter :: kbPerBaryon_to_ergPerKg = 1./ergPerKg_to_kbPerBaryon

   ! integer, save  :: eos_hybVecLen_WL, eos_hybVecLen_HLM

   ! real, dimension(:), allocatable, target, save :: eos_hybData_WL, eos_hybData_HLM

   !! TODO: add declare target directives here?  Not sure if this is necessary here or
   !! if adding target enter/exit data directives after/before allocation/deallocation
   !! is sufficient

end module eos_hybridData
