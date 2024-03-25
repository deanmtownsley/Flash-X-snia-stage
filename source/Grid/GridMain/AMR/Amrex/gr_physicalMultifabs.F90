!!****if* source/Grid/GridMain/AMR/Amrex/gr_physicalMultifabs
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
!!  gr_physicalMultifabs
!!
!! SYNOPSIS
!!
!!  use gr_physicalMultifabs
!!
!! DESCRIPTION 
!!  
!!  Define variables that store physical data as arrays of AMReX
!!  multifabs or as arrays of AMReX flux registers.
!!
!!  In all cases, the array index is the refinement level of the associated
!!  element.  Following AMReX, the level zero corresponds to the coarsest level
!!  of refinement.
!!
!!***

module gr_physicalMultifabs
    use amrex_amr_module,          ONLY : amrex_multifab
    use gr_fluxregister_mod,       ONLY : gr_fluxregister_t

    implicit none

    public :: unk
    public :: facevars

    type(amrex_multifab), allocatable, target :: unk(:)
    type(amrex_multifab), allocatable         :: gr_scratchCtr(:)
    type(amrex_multifab), allocatable         :: facevars(:, :)

    type(amrex_multifab),     allocatable :: fluxes(:, :)
#if N_DIM < 3
    ! Sometimes, pointer arrays need to point somewhere, just anywhere, in order
    ! to have valid code. This is particularly useful when such arrays are
    ! passed as actual arguments to subprograms for corresponding non-pointer
    ! dummy arguments whose values will never be used or set, but which still
    ! need to be present.
    real,save, target :: gr_fakeEmpty4(0,0,0,0) ! fake 4-dimensional array target
#endif
    type(gr_fluxregister_t),  allocatable :: flux_registers(:)

end module gr_physicalMultifabs

