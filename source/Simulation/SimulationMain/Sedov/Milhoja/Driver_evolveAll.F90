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
   use Grid_interface, ONLY: Grid_fillGuardCells
   use Orchestration_interface, ONLY: Orchestration_checkInternalError, &
                                      Orchestration_executeTasks_Cpu
   use Hydro_data, ONLY: hy_gcMaskSize, &
                         hy_gcMask

   !!!!!----- START INSERTION BY CODE GENERATOR
   use dr_hydroAdvance_bundle_mod, ONLY: dr_hydroAdvance_TF_tile_cpu, &
                                         instantiate_hydro_advance_wrapper_C, &
                                         delete_hydro_advance_wrapper_C
   !!!!!----- END INSERTION BY CODE GENERATOR

   implicit none

   character(len=MAX_STRING_LENGTH) :: strBuff(4, 2)
   character(len=15)                :: numToStr
   logical                          :: shortenedDt
   logical                          :: endRun

   !!!!!----- START INSERTION BY CODE GENERATOR
   type(C_PTR)          :: dr_hydroAdvance_wrapper
   integer              :: dr_hydroAdvance_nThreads
   real(MILHOJA_REAL)   :: MH_dt
   integer(MILHOJA_INT) :: MH_ierr

   dr_hydroAdvance_wrapper = C_NULL_PTR 
   !!!!!----- END INSERTION BY CODE GENERATOR

   endRun = .FALSE.

   !!!!!----- START INSERTION BY CODE GENERATOR
   ! RPs are used directly by the Driver and therefore, should be handled at
   ! this level rather than at the level of the code generated for use by
   ! the runtime (i.e., dr_hydroAdvance_bundle_mod).
   CALL RuntimeParameters_get("dr_hydroAdvance_nThreads", &
                              dr_hydroAdvance_nThreads)
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
      MH_dt = REAL(dr_dt, kind=MILHOJA_REAL)
      MH_ierr = instantiate_hydro_advance_wrapper_C(MH_dt, &
                                                    dr_hydroAdvance_wrapper)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
      CALL Orchestration_executeTasks_Cpu(MH_taskFunction=dr_hydroAdvance_TF_tile_cpu, &
                                          prototype_Cptr=dr_hydroAdvance_wrapper, &
                                          nThreads=dr_hydroAdvance_nThreads)
      MH_ierr = delete_hydro_advance_wrapper_C(dr_hydroAdvance_wrapper)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
      dr_hydroAdvance_wrapper = C_NULL_PTR 
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

   CALL Timers_getSummary(MAX(0, dr_nstep - dr_nbegin + 1))
   CALL Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   CALL Logfile_close()
end subroutine Driver_evolveAll

