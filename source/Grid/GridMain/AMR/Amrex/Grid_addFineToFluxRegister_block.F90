!!****if* source/Grid/Grid_addFineToFluxRegister_block
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
!! DESCRIPTION 
!!  Please refer to the documentation in the null implementation of this routine
!!  for more information.
!!
!!***

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif

#include "Simulation.h"
#include "constants.h"

subroutine Grid_addFineToFluxRegister_block(blockDesc, isDensity, coefficient, &
                                           zeroFullRegister)
    use amrex_fort_module,         ONLY : wp => amrex_real
    use amrex_box_module,          ONLY : amrex_box
    use amrex_fab_module,          ONLY : amrex_fab, &
                                          amrex_fab_build, &
                                          amrex_fab_destroy
    use amrex_amrcore_module,      ONLY : amrex_get_finest_level, &
                                          amrex_ref_ratio

    use Driver_interface,          ONLY : Driver_abort
    use Grid_interface,            ONLY : Grid_getGeometry
    use gr_physicalMultifabs,      ONLY : flux_registers
    use Grid_tile,                 ONLY : Grid_tile_t

    implicit none

    type(Grid_tile_t), intent(IN)           :: blockDesc
    logical,           intent(IN), optional :: isDensity(:)
    real,              intent(IN), optional :: coefficient
    logical,           intent(IN), optional :: zeroFullRegister

    integer  :: fine_level
    integer  :: coarse
    integer  :: fine
    integer  :: geometry
    real(wp) :: coef

    real, pointer   :: fluxData(:,:,:,:)
    real, pointer   :: fabData(:,:,:,:)
    type(amrex_box) :: box
    type(amrex_fab) :: fluxFabs(1:NDIM)

    logical :: myZeroFullRegister

    integer :: lo(4)
    integer :: hi(4)

    integer :: axis

    nullify(fluxData)
    nullify(fabData)

    if (NFLUXES < 1) then
        RETURN
    end if

    fine_level = blockDesc%level
    ! FLASH uses 1-based level index / AMReX uses 0-based index
    coarse = fine_level - 2
    fine   = fine_level - 1

    ! The coarsest refinement level is never the fine level of a flux register
    if ((fine <= 0) .OR. (fine > amrex_get_finest_level())) then
#ifdef DEBUG_GRID
        print*,'fine=',fine,', but amrex_get_finest_level() returns',amrex_get_finest_level(),'.'
        call Driver_abort("[Grid_addFineToFluxRegister_block] Invalid level")
#else
        RETURN
#endif
    end if

    if (present(coefficient)) then
        coef = coefficient
    else
        coef = 1.0_wp
    end if

    ! DEV: Determine if this is needed for FLASH5 and implement if necessary. 
    if (present(isDensity)) then
        call Driver_abort("[Grid_addFineToFluxRegister_block] isDensity not implemented")
    end if

    myZeroFullRegister = .FALSE.
    if (present(zeroFullRegister)) then
        myZeroFullRegister = zeroFullRegister
    end if

    call Grid_getGeometry(geometry)

    select case (geometry)
    case (CARTESIAN)
       ! The scaling factor=1/r^(NDIM-1) used here assumes that the refinement
       ! ratio, r, between levels is always 2
       if (amrex_ref_ratio(coarse) /= 2) then
         call Driver_abort("[Grid_addFineToFluxRegister_block] refinement ratio not 2")
       end if

#if   NDIM == 2
       coef = coef * 0.5_wp
#elif NDIM == 3
       coef = coef * 0.25_wp
#endif

       call blockDesc%getDataPtr(fluxData, FLUXX)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       
       box%lo(:) = 1
       box%hi(:) = 1
       box%lo(1:NDIM) = lo(1:NDIM) - 1
       box%hi(1:NDIM) = hi(1:NDIM) - 1
       box%nodal = [.TRUE., .FALSE., .FALSE.]
       call amrex_fab_build(fluxFabs(IAXIS), box, NFLUXES)    
       fabData => fluxFabs(IAXIS)%dataPtr()

       fabData(:, :, :, :) = fluxData(:, :, :, :)
       nullify(fabData)
       call blockDesc%releaseDataPtr(fluxData, FLUXX)

#if   NDIM >= 2
       call blockDesc%getDataPtr(fluxData, FLUXY)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       
       box%lo(:) = 1
       box%hi(:) = 1
       box%lo(1:NDIM) = lo(1:NDIM) - 1
       box%hi(1:NDIM) = hi(1:NDIM) - 1
       box%nodal = [.FALSE., .TRUE., .FALSE.]
       call amrex_fab_build(fluxFabs(JAXIS), box, NFLUXES)    
       fabData => fluxFabs(JAXIS)%dataPtr()

       fabData(:, :, :, :) = fluxData(:, :, :, :)
       nullify(fabData)
       call blockDesc%releaseDataPtr(fluxData, FLUXY)
#endif
#if   NDIM == 3
       call blockDesc%getDataPtr(fluxData, FLUXZ)
       lo(:) = lbound(fluxData)
       hi(:) = ubound(fluxData)
       
       box%lo(:) = 1
       box%hi(:) = 1
       box%lo(1:NDIM) = lo(1:NDIM) - 1
       box%hi(1:NDIM) = hi(1:NDIM) - 1
       box%nodal = [.FALSE., .FALSE., .TRUE.]
       call amrex_fab_build(fluxFabs(KAXIS), box, NFLUXES)    
       fabData => fluxFabs(KAXIS)%dataPtr()

       fabData(:, :, :, :) = fluxData(:, :, :, :)
       nullify(fabData)
       call blockDesc%releaseDataPtr(fluxData, FLUXZ)
#endif

!       call flux_registers(fine)%fineadd(fluxFabs, blockDesc%grid_index, &coef, myZeroFullRegister)

       do axis = 1, NDIM
          call amrex_fab_destroy(fluxFabs(axis))
       end do
    case default
        call Driver_abort("[Grid_addFineToFluxRegister_block] Unsupported geometry")
    end select
end subroutine Grid_addFineToFluxRegister_block

