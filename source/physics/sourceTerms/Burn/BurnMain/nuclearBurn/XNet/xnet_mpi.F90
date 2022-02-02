!***************************************************************************************************
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
! xnet_mpi.f90 06/27/18
! These wrappers are used so that a non-MPI version can coexist with the MPI version.
! Ported from https://github.com/AMReX-Codes/amrex/blob/master/Src/F_BaseLib/parallel.f90 
!***************************************************************************************************

module xnet_mpi
  implicit none
  include 'mpif.h'
  !use mpi
end module xnet_mpi
