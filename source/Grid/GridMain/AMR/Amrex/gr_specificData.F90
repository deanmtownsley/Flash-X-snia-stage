!!****if* source/Grid/GridMain/AMR/Amrex/gr_specificData
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
!!  gr_specificData
!!
!! SYNOPSIS
!!
!!  use gr_specificData
!!
!! DESCRIPTION
!!
!!  Module that contains some private data specific to a Grid implementation.
!!
!!  This includes data that is used by one or more specific implementations,
!!  is probably not of general enough use to be in the Grid_data module,
!!  but is also not part of the internal implementation of an underlying
!!  Grid kernel.
!!
!!  The current Amrex version is derived from the paramesh version, and
!!  should probably change more.
!!***

Module gr_specificData

  implicit none

  integer, save, allocatable :: gr_nToLeft(:)

  ! variables for making visible some Grid metainformation in a PARAMESH-like form
  ! to a friendly IO implementation
  integer, save              :: gr_ioLocalNumBlocks
  integer, save, allocatable :: gr_ioBlkLrefine(:)
  integer, save, allocatable :: gr_ioBlkNodeType(:)
  real,    save, allocatable :: gr_ioBlkCoords(:,:)
  real,    save, allocatable :: gr_ioBlkBsize(:,:)
  real,    save, allocatable :: gr_ioBlkBoundBox(:,:,:)

  logical, save :: gr_enableMaskedGCFill
  integer,save :: gr_globalOffset !stores beginning blk offset for a proc

  logical, save :: gr_bndGCFillNeedsPrimitiveVars

end Module gr_specificData
