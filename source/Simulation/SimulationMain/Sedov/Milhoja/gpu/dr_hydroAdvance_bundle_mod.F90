!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
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
!!
!! @file

#include "constants.h"
#include "Simulation.h"

! The use of this module implies that the runtime is being used.
! Therefore, the Driver unit can assume that precision checks were
! performed by Orchestration_init.  Hence, type casts need not be checked.
#ifndef FLASHX_ORCHESTRATION_MILHOJA
#error "Task function bundles imply the use of the Milhoja runtime"
#endif

!> @details
!! A module that encapsulates all work and specification of data items to which the
!! work should be applied during one invocation of the Orchestration unit.  In
!! other words, the contents of this module specifies implicitly the tasks to be
!! executed by the Orchestration unit.
!!
!! This code should be written entirely by the CODE GENERATOR as part of the
!! full set of bundles needed to execute the timestep recipe.
!!
!! For the purposes of this example, we are imagining that the recipe for our
!! simulation has stipulated that the advancement of the solution due to
!! hydrodynamics effects be computed using both the CPU and GPU simultaneously.  For
!! example, these two routines could be mapped onto the CPU/GPU data parallel
!! thread team configuration.  In addition, we are imagining that the CODE
!! GENERATOR has been instructed to use OpenAcc directives to offload
!! computation to the GPU.
!!
!! Note that this module does not specify how the contents are mapped onto 
!! a thread team configuration.
!!
!! @todo Add in all metadata related to the file, the code generator, and the
!!       Simulation that this was created for.  This could include, for example,
!!         - timestamp of creation
!!         - git commit and state of repo at time of generation
!!         - state of files in repo that provide information used to generate
!!           this file and diffs if altered.  Do such files exist?
!!         - machine on which run and username
!!         - name/version of code generator
!!         - setup call
!!         - name/version of recipe that these are to be used in
!! @todo Is this a good design?  Main driver is likely to keep CODE GENERATORS
!!       simple and to have products to generally-useful.  For instance, we
!!       don't want the generated code to be too tightly-linked to Flash-X.
module dr_hydroAdvance_bundle_mod
    implicit none
    private

#ifdef ORCHESTRATION_USE_GPUS
    public :: dr_hydroAdvance_packet_gpu_oacc
    public :: dr_hydro_advance_packet_oacc_tf
#endif

   interface
      subroutine dr_hydro_advance_packet_oacc_tf(C_tId, C_dataItemPtr) bind(c)
         use iso_c_binding,     ONLY : C_PTR
         use milhoja_types_mod, ONLY : MILHOJA_INT
         integer(MILHOJA_INT), intent(IN), value :: C_tId
         type(C_PTR),          intent(IN), value :: C_dataItemPtr
      end subroutine dr_hydro_advance_packet_oacc_tf
   end interface

contains

#if defined(ORCHESTRATION_USE_GPUS) && defined(ORCHESTRATION_OPENACC_OFFLOAD)
    !> @brief GPU-based variant of the hydro advance task function
    !!
    !! @details
    !! This routine is *not* a task function that is called directly by a pthread
    !! executing C++ code.  Rather the task function is C++ code
    !! (dr_hydro_advance_packet_oacc_tf) that knows how to pick out the contents of
    !! the given data packet and map it indirectly onto this function's argument
    !! list via a C-to-Fortran reinterpreation layer
    !! (dr_hydro_advance_packet_oacc_C2F).  Despite this, this routine is still
    !! part of the Task Function/Data Packet Fortran/C++ interoperability design
    !! and should satisy that design's requirements.  This justifies the fact that
    !! the interface exposes iso_c_binding types.
    !!
    !! This code uses OpenACC asynchronous queues so that code is only run once the
    !! given data packet has arrived in the GPU's memory system.  It is assumed
    !! that the data packet transfer was scheduled using dataQ_h.
    !!
    !! The _h variable suffix indicates that the variable is located in the
    !! host's memory system; the _d suffix, in the GPU's memory system. The
    !! _all_d variable suffix indicates that the variable is an array of objects
    !! (themselves potentially arrays) collected across all data items in the
    !! given packet.
    !!
    !! @todo The lo_loc_n/hi_loc_n local indices are presently used as a workaround
    !! associated with a Summit/NVHPC compiler error that is not yet understood.
    !! It arises when calling the *_gpu_oacc routines called here.  It would be
    !! best to get rid of these and just pass in lo/hi global indices directly.
    !! Else, the runtime's latency hiding is potentially lost due to the memory
    !! allocation.
    !! @todo The 2D and 3D kernel launching scheme has not been optimized for
    !! Fortran/Flash-X.  Rather, these were ported directly from C++ code.
    !! Study the best scenario for each setup.
    !! @todo Fuse updateSolution & EoS?
    subroutine dr_hydroAdvance_packet_gpu_oacc(C_packet_h,         &
                                               dataQ_h,            &
                                               queue2_h, queue3_h, &
                                               nTiles_d, dt_d,     &
                                               deltas_all_d,       &
                                               lo_all_d, hi_all_d, &
                                               loGC_all_d,         &
                                               U_all_d,            &
                                               auxC_all_d,         &
                                               faceX_all_d,        &
                                               faceY_all_d,        &
                                               faceZ_all_d)
        use iso_c_binding, ONLY : C_PTR
        use openacc

        use milhoja_types_mod,               ONLY : MILHOJA_INT
        use Orchestration_Interface,         ONLY : Orchestration_checkInternalError
        use DataPacket_gpu_tf_hydro_C2F_mod, ONLY : release_hydro_advance_extra_queue_C 
        use dr_cg_hydroAdvance_mod,          ONLY : hy_computeSoundSpeedHll_gpu_oacc, &
                                                    hy_computeFluxesHll_X_gpu_oacc,   &
                                                    hy_computeFluxesHll_Y_gpu_oacc,   &
                                                    hy_computeFluxesHll_Z_gpu_oacc,   &
                                                    hy_updateSolutionHll_gpu_oacc,    &
                                                    eos_idealGammaDensIe_gpu_oacc

        implicit none

        !$acc routine (hy_computeSoundSpeedHll_gpu_oacc) vector
        !$acc routine (hy_computeFluxesHll_X_gpu_oacc)   vector
        !$acc routine (hy_computeFluxesHll_Y_gpu_oacc)   vector
        !$acc routine (hy_computeFluxesHll_Z_gpu_oacc)   vector
        !$acc routine (hy_updateSolutionHll_gpu_oacc)    vector
        !$acc routine (eos_idealGammaDensIe_gpu_oacc)    vector

        type(C_PTR),                   intent(IN)    :: C_packet_h
        integer(kind=acc_handle_kind), intent(IN)    :: dataQ_h
        integer(kind=acc_handle_kind), intent(IN)    :: queue2_h
        integer(kind=acc_handle_kind), intent(IN)    :: queue3_h
        integer,                       intent(IN)    :: nTiles_d
        real,                          intent(IN)    :: dt_d
        real,                          intent(IN)    :: deltas_all_d(:, :)
        integer,                       intent(IN)    :: lo_all_d(:, :)
        integer,                       intent(IN)    :: hi_all_d(:, :)
        integer,                       intent(IN)    :: loGC_all_d(:, :)
        real,                          intent(INOUT) :: U_all_d(:, :, :, :, :)
        real,                          intent(OUT)   :: auxC_all_d(:, :, :, :)
        real,                          intent(OUT)   :: faceX_all_d(:, :, :, :, :)
        real,                          intent(OUT)   :: faceY_all_d(:, :, :, :, :)
        real,                          intent(OUT)   :: faceZ_all_d(:, :, :, :, :)

        integer :: lo_loc_n(1:MDIM)
        integer :: hi_loc_n(1:MDIM)
        integer :: n

        integer(MILHOJA_INT) :: MH_id
        integer(MILHOJA_INT) :: MH_ierr

        !$acc  data create(lo_loc_n, hi_loc_n)                       &
        !$acc&      deviceptr(nTiles_d, dt_d, deltas_all_d,          &
        !$acc&                lo_all_d, hi_all_d, loGC_all_d,        &
        !$acc&                faceX_all_d, faceY_all_d, faceZ_all_d, &
        !$acc&                auxC_all_d, U_all_d)

        !----- ADVANCE SOLUTION
        ! Update UNK data on interiors only
        !   * It is assumed that the GC are filled already
        !   * No tiling for now means that computing fluxes and updating the
        !     solution can be fused and the full advance run independently on each
        !     block and in place.

        !----- COMPUTE FLUXES
        ! The X, Y, and Z fluxes each depend on the speed of sound, but can
        ! be computed independently and therefore concurrently.
        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(dataQ_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            CALL hy_computeSoundSpeedHll_gpu_oacc(lo_loc_n, hi_loc_n,     &
                                                  U_all_d(:, :, :, :, n), &
                                                  auxC_all_d(:, :, :, n))
        end do
        !$acc end parallel loop

#if   NDIM == 1
        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(dataQ_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            CALL hy_computeFluxesHll_X_gpu_oacc(dt_d,                       &
                                                lo_loc_n, hi_loc_n,         &
                                                deltas_all_d(:, n),         &
                                                U_all_d(:, :, :, :, n),     &
                                                auxC_all_d(:, :, :, n),     &
                                                faceX_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        ! No need for barrier since only one kernel launched in 1D case
#elif NDIM == 2
        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(dataQ_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            ! It seems like for small 2D blocks, fusing kernels is more
            ! efficient than fusing actions (i.e. running the two kernels
            ! concurrently).  Too much work for the GPU?  Too much overhead
            ! from the stream sync (i.e. OpenACC wait)?
            CALL hy_computeFluxesHll_X_gpu_oacc(dt_d,                       &
                                                lo_loc_n, hi_loc_n,         &
                                                deltas_all_d(:, n),         &
                                                U_all_d(:, :, :, :, n),     &
                                                auxC_all_d(:, :, :, n),     &
                                                faceX_all_d(:, :, :, :, n))
            CALL hy_computeFluxesHll_Y_gpu_oacc(dt_d,                       &
                                                lo_loc_n, hi_loc_n,         &
                                                deltas_all_d(:, n),         &
                                                U_all_d(:, :, :, :, n),     &
                                                auxC_all_d(:, :, :, n),     &
                                                faceY_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        ! No need for barrier since all kernels are launched on the same
        ! queue for 2D case.
#elif NDIM == 3
        ! Wait for computation of speed of sound first and then launch
        ! the flux computations concurrently.
        !$acc wait(dataQ_h)

        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(queue2_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            CALL hy_computeFluxesHll_X_gpu_oacc(dt_d,                       &
                                                lo_loc_n, hi_loc_n,         &
                                                deltas_all_d(:, n),         &
                                                U_all_d(:, :, :, :, n),     &
                                                auxC_all_d(:, :, :, n),     &
                                                faceX_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(queue3_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            CALL hy_computeFluxesHll_Y_gpu_oacc(dt_d,                       &
                                                lo_loc_n, hi_loc_n,         &
                                                deltas_all_d(:, n),         &
                                                U_all_d(:, :, :, :, n),     &
                                                auxC_all_d(:, :, :, n),     &
                                                faceY_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(dataQ_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            CALL hy_computeFluxesHll_Z_gpu_oacc(dt_d,                       &
                                                lo_loc_n, hi_loc_n,         &
                                                deltas_all_d(:, n),         &
                                                U_all_d(:, :, :, :, n),     &
                                                auxC_all_d(:, :, :, n),     &
                                                faceZ_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        !$acc wait(queue2_h, queue3_h)
        MH_id = INT(2, kind=MILHOJA_INT)
        MH_ierr = release_hydro_advance_extra_queue_C(C_packet_h, MH_id)
        CALL Orchestration_checkInternalError("dr_hydroAdvance_packet_gpu_oacc", MH_ierr)
        MH_id = INT(3, kind=MILHOJA_INT)
        MH_ierr = release_hydro_advance_extra_queue_C(C_packet_h, MH_id)
        CALL Orchestration_checkInternalError("dr_hydroAdvance_packet_gpu_oacc", MH_ierr)
#endif

        !----- UPDATE SOLUTIONS IN PLACE
        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(dataQ_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            ! NOTE: If NDIM < 3, then some of the face[YZ]_all_d will be garbage.
            !       We therefore assume that this routine will not use
            !       those fluxes associated with axes "above" NDIM.  We still
            !       want indexing into it by tile to yield a valid array.
            CALL hy_updateSolutionHll_gpu_oacc(lo_loc_n, hi_loc_n,         &
                                               faceX_all_d(:, :, :, :, n), &
                                               faceY_all_d(:, :, :, :, n), &
                                               faceZ_all_d(:, :, :, :, n), &
                                               U_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        !$acc  parallel loop gang default(none)               &
        !$acc&                    private(lo_loc_n, hi_loc_n) &
        !$acc&                    async(dataQ_h)
        do n = 1, nTiles_d
            lo_loc_n(:) = lo_all_d(:, n) - loGC_all_d(:, n) + 1
            hi_loc_n(:) = hi_all_d(:, n) - loGC_all_d(:, n) + 1

            CALL eos_idealGammaDensIe_gpu_oacc(lo_loc_n, hi_loc_n,     &
                                               U_all_d(:, :, :, :, n))
        end do
        !$acc end parallel loop

        !$acc wait(dataQ_h)

        !$acc end data
    end subroutine dr_hydroAdvance_packet_gpu_oacc
#endif

end module dr_hydroAdvance_bundle_mod

