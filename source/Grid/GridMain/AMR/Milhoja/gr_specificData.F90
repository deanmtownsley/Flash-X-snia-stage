!> @file
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock

!> @ingroup GridMilhoja
!!
!! @brief
!! Milhoja-specific local data variables
!!
!! @details
!! Module that contains some private data specific to a Grid implementation that is
!! probably not of general enough use to be in the Grid_data module, and is also
!! not part of the internal implementation of an underlying Grid kernel.
!!
!! @todo This was seeded from the AMReX version.  Need to determine if we need
!! to build up anything beyond this.  So far, I have confirmed that the current
!! contents are necessary.
Module gr_specificData
  implicit none

  ! These are used to communicate to the IO unit data needed for checkpointing.
  ! Therefore the IO unit is a "friend" of Grid in the sense that it can access
  ! these variables.  It is not intended that these be used by any other code
  ! and they should only be used by IO if Grid_sendOutputData has been called.
  !
  ! variables for making visible some Grid metainformation in a PARAMESH-like form
  ! to a friendly IO implementation
  integer, save              :: gr_ioLocalNumBlocks
  integer, save              :: gr_ioGlobalNumBlocks
  integer, save, allocatable :: gr_ioBlkLrefine(:)
  integer, save, allocatable :: gr_ioBlkNodeType(:)
  real,    save, allocatable :: gr_ioBlkCoords(:, :)
  real,    save, allocatable :: gr_ioBlkBsize(:, :)
  real,    save, allocatable :: gr_ioBlkBoundBox(:, :, :)
end Module gr_specificData

