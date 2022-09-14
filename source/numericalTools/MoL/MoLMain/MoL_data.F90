!!****if* source/numericalTools/MoL/MoLMain/MoL_data
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
!!  NAME
!!    MoL_data
!!
!!  SYNOPSIS
!!    use MoL_data
!!
!!  DESCRIPTION
!!    Stores data for MoL
!!
!!  NOTES
!!
!!***
module MoL_data

   implicit none

   ! ml_nscratch is the required number of memory-levels (e.g. intermediate
   ! stages) required by an implementation
   ! ml_nscratch_total adds additional levels for default memory-levels provided
   ! by the base implementation (e.g. RHS and INITIAL)
   integer, save :: ml_nscratch, ml_nscratch_total

   ! Current MPI rank (mainly used for messaging)
   integer, save :: ml_mpiRank

   ! Verbosity for messaging
   integer, save :: ml_verbosity
   logical, save :: ml_abortOnWarn

end module MoL_data
