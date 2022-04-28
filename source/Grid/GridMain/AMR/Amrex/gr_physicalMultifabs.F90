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
    public :: facevarx
    public :: facevary
    public :: facevarz

    type(amrex_multifab), allocatable, target :: unk(:)
    type(amrex_multifab), allocatable         :: facevarx(:)
    type(amrex_multifab), allocatable         :: facevary(:)
    type(amrex_multifab), allocatable         :: facevarz(:)
    type(amrex_multifab), allocatable         :: gr_scratchCtr(:)

    type(amrex_multifab),     allocatable :: fluxes(:, :)
    type(gr_fluxregister_t),  allocatable :: flux_registers(:)

end module gr_physicalMultifabs

