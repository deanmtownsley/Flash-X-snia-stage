#include "constants.h"
#include "Simulation.h"

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
!! Wrapper functions that allow Simulation_initBlock code to be run by Milhoja.
!! In the future, this module *might* include interfaces for using the runtime
!! with different thread team configurations (TBD).
!! 
!! All routines in this module require that Simulation_initBlock developers write
!! initial condition data to the interiors and GCs of all given tiles.  After
!! calling Simulation_initBlock on the given tile, each routine writes EoS initial
!! conditions to the interior and GCs of the block.  The ICs must be written to the
!! GCs since this routine must leave the blocks ready for immediate use by
!! refinement/derefinement routines and these might choose, for example, to refine
!! a block based solely on the presence of a fine feature in the GCs.  Therefore,
!! at the very least, initial condition data must be written for all cell-centered
!! variables that serve as refinement variables or that are needed to compute EoS
!! variables that are refinement variables.
!!
!! @todo Can we have GPU-only and CPU/GPU data parallel versions?
!! @todo Flash-X does not require that Simulation_initBlock write GC data.
!!       Therefore, this needs to conform to that non-constraint.
!! @todo Would it be better to allow each simulation to specify if
!!       Simulation_initBlock should set ICs in the GCs?
module gr_initDomain_mod
    implicit none
    private

    public :: gr_initBlock_tile_cpu
    public :: gr_initBlock_wrapper_cpu
    public :: gr_instantiate_wrapper_C
    public :: gr_delete_wrapper_C

    !!!!!----- INTERFACES TO C-LINKAGE C++ FUNCTIONS
    ! The C-to-Fortran interoperability layer
    interface
        function gr_instantiate_wrapper_C(C_wrapper) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            type(C_PTR),         intent(IN) :: C_wrapper
            integer(MILHOJA_INT)            :: C_ierr
        end function gr_instantiate_wrapper_C

        function gr_delete_wrapper_C(C_wrapper) result(C_ierr) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            type(C_PTR),         intent(IN), value :: C_wrapper
            integer(MILHOJA_INT)                   :: C_ierr
        end function gr_delete_wrapper_C
    end interface

contains

    !> Grid_initDomain should use this function to compute the ICs for a single,
    !! given tile using the CPU either with or without the Milhoja runtime.
    !!
    !! This wrapper patches the contents of the given data item onto the
    !! argument list of Simulation_initBlock.
    !!
    !! Since this can be registered with the runtime as a task function, it is
    !! part of the Fortran/C++ interoperability layer.  As a result it is
    !! acceptable that it have a C-linkage-compatible interface.
    !!
    !! @todo Look at AMReX implementation of gr_initNewLevelCallback and
    !!       implement all actions needed there for setting ICs.
    !!
    !! @param C_dataItemPtr   A C pointer to the tile whose data should be
    !!                        set to the initial conditions using
    !!                        Simulation_initBlock
    subroutine gr_initBlock_tile_cpu(C_dataItemPtr) bind(c)
        use iso_c_binding,        ONLY : C_PTR

        use milhoja_types_mod,    ONLY : MILHOJA_INT

        use Grid_data,            ONLY : gr_eosModeInit
        use Grid_tile,            ONLY : Grid_tile_t, &
                                         Grid_tile_fromMilhojaTilePtr
        use Eos_interface,        ONLY : Eos_multiDim
        use Simulation_interface, ONLY : Simulation_initBlock

        type(C_PTR),          intent(IN), value :: C_dataItemPtr

        real, contiguous, pointer :: initData(:,:,:,:)

        type(Grid_tile_t) :: tileDesc
        integer           :: i, j, k
        integer           :: var

        NULLIFY(initData)

        CALL Grid_tile_fromMilhojaTilePtr(C_dataItemPtr, tileDesc)

        ! We need to zero data in case we reuse blocks from previous levels
        ! but don't initialize all data in Simulation_initBlock... in particular
        ! the total vs. internal energies can cause problems in the eos call that 
        ! follows.  This includes zeroing the data in the grown tile's
        ! guardcells.
        associate(loGC => tileDesc%grownLimits(LOW,  :), &
                  hiGC => tileDesc%grownLimits(HIGH, :))
            CALL tileDesc%getDataPtr(initData, CENTER)
            do           var = UNK_VARS_BEGIN, UNK_VARS_END
                do         k = loGC(KAXIS), hiGC(KAXIS)
                    do     j = loGC(JAXIS), hiGC(JAXIS)
                        do i = loGC(IAXIS), hiGC(IAXIS)
                            initData(i, j, k, var) = 0.0
                        end do
                    end do
                end do
            end do

            CALL Simulation_initBlock(initData, tileDesc)
            CALL Eos_multiDim(gr_eosModeInit, tileDesc%grownLimits, initData)

            CALL tileDesc%releaseDataPtr(initData, CENTER)
        end associate
    end subroutine gr_initBlock_tile_cpu

    subroutine gr_initBlock_wrapper_cpu(C_threadId, C_dataItemPtr) bind(c)
        use iso_c_binding,        ONLY : C_PTR, &
                                         C_NULL_PTR

        use milhoja_types_mod,   ONLY : MILHOJA_INT
        use milhoja_tile_mod,    ONLY : milhoja_tile_from_wrapper_C
        use gr_milhojaInterface, ONLY : gr_checkMilhojaError

        integer(MILHOJA_INT), intent(IN), value :: C_threadId
        type(C_PTR),          intent(IN), value :: C_dataItemPtr

        integer(MILHOJA_INT) :: MH_ierr

        type(C_PTR) :: C_tilePtr

        C_tilePtr = C_NULL_PTR

        MH_ierr = milhoja_tile_from_wrapper_C(C_dataItemPtr, C_tilePtr)
        CALL gr_checkMilhojaError("gr_initBlock_wrapper_cpu", MH_ierr)

        CALL gr_initBlock_tile_cpu(C_tilePtr)
    end subroutine gr_initBlock_wrapper_cpu

end module gr_initDomain_mod
