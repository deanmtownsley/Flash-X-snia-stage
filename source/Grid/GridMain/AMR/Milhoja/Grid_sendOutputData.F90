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

#include "constants.h"
#include "Simulation.h"

!> @ingroup GridMilhoja
!! @stubref{Grid_sendOutputData}
!!
!! @brief Concrete implementation of Grid_sendOutputData
!!
!! @attention
!! This implementation was built up quick & dirty from the AMReX routine.
!! In particular, the current form is only sufficient for pseudo-UG.
!!
!! @todo Complete full implementation
subroutine Grid_sendOutputData()
    use Grid_data,           ONLY : gr_str_geometry, &
                                    gr_meshComm, &
                                    gr_meshNumProcs
    use gr_specificData,     ONLY : gr_ioLocalNumBlocks, &
                                    gr_ioGlobalNumBlocks
    use gr_milhojaInterface, ONLY : gr_updateDataForIo
    use IO_interface,        ONLY : IO_setScalar

    implicit none

    include "Flashx_mpi.h"

    ! Put NXB, NYB and NZB into a saved variable to prevent problems with
    ! an xlf "feature."
    integer, parameter :: LOCAL_NXB      = NXB
    integer, parameter :: LOCAL_NYB      = NYB
    integer, parameter :: LOCAL_NZB      = NZB
    integer, parameter :: DIMENSIONALITY = NDIM

    integer :: localNumBlocks
    integer :: nToLeft(0:gr_meshNumProcs - 1)
    integer :: ierr
    integer :: i

    !set the scalars for the grid unit
    CALL IO_setScalar("nxb", LOCAL_NXB)
    CALL IO_setScalar("nyb", LOCAL_NYB)
    CALL IO_setScalar("nzb", LOCAL_NZB)

    CALL IO_setScalar("geometry", gr_str_geometry)
    CALL IO_setScalar("dimensionality", DIMENSIONALITY)

    ! This (re)computes gr_ioLocalNumBlocks and gr_ioBlkBLAH arrays.
    CALL gr_updateDataForIo()

    ! Get the local number of blocks from everybody
    ! to find total number of blocks
    localNumBlocks = gr_ioLocalNumBlocks
    CALL MPI_Allgather(localNumBlocks, 1, MPI_INTEGER, nToLeft, &
                       1, MPI_INTEGER, gr_meshComm, ierr)

    gr_ioGlobalNumBlocks = 0
    do i = 0, gr_meshNumProcs - 1
        gr_ioGlobalNumBlocks = gr_ioGlobalNumBlocks + nToLeft(i)
    end do

    CALL IO_setScalar("globalNumBlocks", gr_ioGlobalNumBlocks)
end subroutine Grid_sendOutputData

