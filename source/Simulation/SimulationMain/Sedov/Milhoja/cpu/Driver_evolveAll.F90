!> @copyright Copyright 2024 UChicago Argonne, LLC and contributors
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

!> @details
!! This is a barebones version of a Driver_evolveAll that is needed to stand-up an
!! official test of the Milhoja Orchestration unit.  One goal is for this to be
!! implemented in such a way that using this to get Fortran-based performance
!! results will yield results that can be reasonably compared against similar
!! C++-based performance results in a publication.
!!
!! @note Names should be chosen so that recipe creators can read this file easily.
!! I don't know if our goal should be for recipe creators to be able to manually
!! alter this file easily.  Ideally, they should be able to just read this file,
!! adjust the recipe accordingly, regenerate with the offline toolchain, and diff
!! the new/old versions if desired.
!! @note I like the idea of using the dummy argument names in the call to the
!! Orchestration unit to make explicit here the mapping.  Detailed names for the
!! actual arguments will hopefully avoid having to open up many files at once just
!! to digest the Fortran expression here of the timestep recipe.
!! @note Johann to determine how to write this to file.  Is using a module
!! for each bundle easiest/cleanest?  If so, is it OK to violate standard Flash-X
!! file design rules?  Does it make sense that this data is owned by Driver since
!! any bundle might include task functions from different physics units?
!! @note The design of the use of runtime parameters for configuration of invocations
!! of the Orchestration unit shall not preclude the possibility of dynamically
!! altering these values from one timestep to the next using acquired performance
!! results and predictive models.
!!
!! @todo Remove from repo once offline toolchain is in place.
!! @todo Compare against C++ version before beginning study.
subroutine Driver_evolveAll()
   use iso_c_binding, ONLY : C_PTR, &
                             C_NULL_PTR
 
   use milhoja_types_mod,   ONLY : MILHOJA_INT, &
                                   MILHOJA_REAL

   use Driver_data, ONLY: dr_nstep, &
                          dr_nbegin, &
                          dr_nend, &
                          dr_simTime, &
                          dr_dt, &
                          dr_dtOld, &
                          dr_dtNew, &
                          dr_tmax, &
                          dr_globalMe
   use Driver_interface, ONLY: Driver_computeDt
   use Logfile_interface, ONLY: Logfile_stamp, &
                                Logfile_close
   use Timers_interface, ONLY: Timers_start, &
                               Timers_stop, &
                               Timers_getSummary
   use RuntimeParameters_interface, ONLY: RuntimeParameters_get
   use IO_interface, ONLY: IO_output, &
                           IO_outputFinal
   use Grid_interface, ONLY: Grid_fillGuardCells, &
                             Grid_getTileIterator, &
                             Grid_releaseTileIterator
   use Grid_iterator, ONLY : Grid_iterator_t
   use Grid_tile,     ONLY : Grid_tile_t
   use Orchestration_interfaceTypeDecl, ONLY: Orchestration_tileCInfo_t
   use Orchestration_interface, ONLY: Orchestration_checkInternalError
#ifdef RUNTIME_USES_TILEITER
   use Orchestration_interface, ONLY: Orchestration_executeTasks_Cpu
#else
   use Orchestration_interface, ONLY: Orchestration_setupPipelineForCpuTasks, &
                                      Orchestration_pushTileToPipeline, &
                                      Orchestration_teardownPipeline
#endif
   use Hydro_data, ONLY: hy_gcMaskSize, &
                         hy_gcMask, &
                         hy_eosModeAfter

   !!!!!----- START INSERTION BY CODE GENERATOR
   use cpu_tf_hydro_mod,          ONLY : cpu_tf_hydro_Cpp2C
   use Tile_cpu_tf_hydro_C2F_mod, ONLY : new_hydro_advance_wrapper_C, &
                                         delete_hydro_advance_wrapper_C, &
                                         acquire_scratch_wrapper_C, &
                                         release_scratch_wrapper_C
   !!!!!----- END INSERTION BY CODE GENERATOR

   implicit none

   character(len=MAX_STRING_LENGTH) :: strBuff(4, 2)
   character(len=15)                :: numToStr
   logical                          :: shortenedDt
   logical                          :: endRun
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t)     :: tileDesc
   type(Orchestration_tileCInfo_t) :: cInfo

   !!!!!----- START INSERTION BY CODE GENERATOR
   type(C_PTR)          :: cpu_tf_hydro_wrapper
   integer              :: cpu_tf_hydro_nThreads
   real(MILHOJA_REAL)   :: MH_dt
   integer(MILHOJA_INT) :: MH_eosMode
   integer(MILHOJA_INT) :: MH_ierr

   cpu_tf_hydro_wrapper = C_NULL_PTR 

   ! Acquire persistent milhoja scratch
   MH_ierr = acquire_scratch_wrapper_C()
   CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
   !!!!!----- END INSERTION BY CODE GENERATOR

   endRun = .FALSE.

   !!!!!----- START INSERTION BY CODE GENERATOR
   ! RPs are used directly by the Driver and therefore should be handled at
   ! this level rather than at the level of the code generated for use by
   ! the runtime (i.e., dr_hydroAdvance_bundle_mod).
   CALL RuntimeParameters_get("cpu_tf_hydro_nThreads", &
                              cpu_tf_hydro_nThreads)
   !!!!!----- END INSERTION BY CODE GENERATOR

   CALL Logfile_stamp('Entering evolution loop', '[Driver_evolveAll]')
   CALL Timers_start("evolution")

   do dr_nstep = dr_nBegin, dr_nend
      CALL dr_shortenLastDt(dr_dt, dr_simTime, dr_tmax, shortenedDt, 1)
      if (dr_globalMe == MASTER_PE) then
         write (numToStr(1:), '(I10)') dr_nstep
         write (strBuff(1, 1), "(A)") "n"
         write (strBuff(1, 2), "(A)") trim(adjustl(numToStr))

         write (numToStr(1:), "(1PE12.6)") dr_simTime
         write (strBuff(2, 1), "(A)") "t"
         write (strBuff(2, 2), "(A)") trim(adjustl(numToStr))

         write (numToStr(1:), "(1PE12.6)") dr_dt
         write (strBuff(3, 1), "(A)") "dt"
         write (strBuff(3, 2), "(A)") trim(adjustl(numToStr))

         CALL Logfile_stamp(strBuff(1:3, :), 3, 2, "step")
      end if
      dr_simTime = dr_simTime + dr_dt

      !!!!!----- START INSERTION BY CODE GENERATOR
      ! This is the timestep for now as dictated by the recipe and rendered as
      ! Fortran code by the CODE GENERATOR.
      !
      ! This maps the tasks defined implicitly by the bundle onto the chosen
      ! thread team configuration.  In this case, the CPU-only config.  The
      ! mapping should have been specified originally in the recipe.
      CALL Grid_fillGuardCells(CENTER, ALLDIR, &
                               doEos=.TRUE., &
                               maskSize=hy_gcMaskSize, &
                               mask=hy_gcMask, &
                               makeMaskConsistent=.TRUE., &
                               selectBlockType=LEAF, &
                               doLogMask=.TRUE.)
      ! We need to instantiate a new TileWrapper prototype at each step since
      ! dt can change from one step to the next.
      MH_dt = REAL(dr_dt, kind=MILHOJA_REAL)
      MH_eosMode = INT(hy_eosModeAfter, kind=MILHOJA_INT)
      MH_ierr = new_hydro_advance_wrapper_C(MH_dt, MH_eosMode, &
                                            cpu_tf_hydro_wrapper)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
#ifdef RUNTIME_USES_TILEITER
      CALL Orchestration_executeTasks_Cpu(MH_taskFunction=cpu_tf_hydro_Cpp2C, &
                                          prototype_Cptr=cpu_tf_hydro_wrapper, &
                                          nThreads=cpu_tf_hydro_nThreads)
#else
      CALL Orchestration_setupPipelineForCpuTasks(MH_taskFunction=cpu_tf_hydro_Cpp2C, &
                                          nThreads=cpu_tf_hydro_nThreads)

      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while(itor%isValid())
         call itor%currentTile(tileDesc)
         call tileDesc%fillTileCInfo(cInfo)
         
         CALL Orchestration_pushTileToPipeline(prototype_Cptr=cpu_tf_hydro_wrapper, &
                                          nThreads=cpu_tf_hydro_nThreads, &
                                          tileCInfo=cInfo)

         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)

      CALL Orchestration_teardownPipeline(nThreads=cpu_tf_hydro_nThreads)
#endif
      MH_ierr = delete_hydro_advance_wrapper_C(cpu_tf_hydro_wrapper)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
      cpu_tf_hydro_wrapper = C_NULL_PTR 
      !!!!!----- END INSERTION BY CODE GENERATOR

      dr_dtOld = dr_dt
      call Driver_computeDt(dr_nbegin, dr_nstep, &
                            dr_simTime, dr_dtOld, dr_dtNew)
      dr_dt = dr_dtNew

      call IO_output(dr_simTime, dr_dt, dr_nstep + 1, dr_nbegin, endRun, &
                     CHECKPOINT_FILE_ONLY)

      if (shortenedDt) then
         exit
      end if

      if (dr_simTime >= dr_tmax) then
         if (dr_globalMe == MASTER_PE) then
            print *, "exiting: reached max SimTime"
         end if
         exit
      end if
   end do
   dr_nstep = min(dr_nstep, dr_nend)

   CALL Timers_stop("evolution")
   CALL Logfile_stamp('Exiting evolution loop', '[Driver_evolveAll]')

   call IO_outputFinal()

   !!!!!----- START INSERTION BY CODE GENERATOR
   ! Release persistent milhoja scratch
   MH_ierr = release_scratch_wrapper_C()
   CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
   !!!!!----- END INSERTION BY CODE GENERATOR

   CALL Timers_getSummary(MAX(0, dr_nstep - dr_nbegin + 1))
   CALL Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   CALL Logfile_close()
end subroutine Driver_evolveAll

