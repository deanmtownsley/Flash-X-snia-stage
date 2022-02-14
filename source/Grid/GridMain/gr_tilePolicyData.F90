!! DEV: This module is currently not useful. - KW
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
!! ... except for (sometimes) producing useful console messages about tile sizes.
#include "Simulation.h"
#include "constants.h"

module gr_tilePolicyData
  implicit none
  integer, save :: gr_tileNxt
#if NDIM > 1
  integer, save :: gr_tileNyt
#else
  integer, parameter :: gr_tileNyt = 1
#endif
#if NDIM > 2
  integer, save :: gr_tileNzt
#else
  integer, parameter :: gr_tileNzt = 1
#endif

  integer,save :: gr_tileNumTiles

contains
  subroutine gr_initTilePolicy
!!$    use RuntimeParameters_interface, ONLY: RuntimeParameters_get
    use Grid_data, ONLY: gr_enableTiling
    use Grid_data, ONLY: gr_tileSize
    use Grid_data, ONLY: gr_meshMe
    implicit none
!!$    call RuntimeParameters_get("gr_tileNxt", gr_tileNxt)
    gr_tileNxt = max(1,(NXB+gr_tileSize(IAXIS)-1)/gr_tileSize(IAXIS))
#if NDIM > 1
!!$    call RuntimeParameters_get("gr_tileNyt", gr_tileNyt)
    gr_tileNyt = max(1,(NYB+gr_tileSize(JAXIS)-1)/gr_tileSize(JAXIS))
#endif
#if NDIM > 2
!!$    call RuntimeParameters_get("gr_tileNzt", gr_tileNzt)
    gr_tileNzt = max(1,(NZB+gr_tileSize(KAXIS)-1)/gr_tileSize(KAXIS))
#endif

    gr_tileNumTiles = gr_tileNxt
#if NDIM > 1
    gr_tileNumTiles = gr_tileNumTiles * gr_tileNyt
#endif
#if NDIM > 2
    gr_tileNumTiles = gr_tileNumTiles * gr_tileNzt
#endif

#ifdef FIXEDBLOCKSIZE
    if (gr_enableTiling) then
       if (gr_meshMe==MASTER_PE) print*,'Number of tiles per block in X-direction:', gr_tileNxt
#if NDIM > 1
       if (gr_meshMe==MASTER_PE) print*,'Number of tiles per block in Y-direction:', gr_tileNyt
#endif
#if NDIM > 2
       if (gr_meshMe==MASTER_PE) print*,'Number of tiles per block in Z-direction:', gr_tileNzt
#endif
    end if
#endif

  end subroutine gr_initTilePolicy
end module gr_tilePolicyData
