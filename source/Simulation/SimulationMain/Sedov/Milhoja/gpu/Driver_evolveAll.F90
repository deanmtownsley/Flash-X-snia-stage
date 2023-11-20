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
!! @todo Compare against C++ version before beginning study.
!! @todo Remove runtime reset once memory manager in place.
subroutine Driver_evolveAll()
   use iso_c_binding, ONLY : C_PTR, &
                             C_NULL_PTR

   use milhoja_types_mod,   ONLY : MILHOJA_INT, &
                                   MILHOJA_REAL
   use milhoja_runtime_mod, ONLY : milhoja_runtime_reset

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
   use Orchestration_interface, ONLY: Orchestration_executeTasks_Gpu, &
                                      Orchestration_checkInternalError
   use Hydro_data, ONLY: hy_gcMaskSize, &
                         hy_gcMask

   !!!!!----- START INSERTION BY CODE GENERATOR
   use DataPacket_gpu_tf_hydro_C2F_mod, ONLY : instantiate_hydro_advance_packet_C, &
                                               delete_hydro_advance_packet_C
   use gpu_tf_hydro_mod,                ONLY : gpu_tf_hydro
   !!!!!----- END INSERTION BY CODE GENERATOR

   implicit none

   character(len=MAX_STRING_LENGTH) :: strBuff(4, 2)
   character(len=15)                :: numToStr
   logical                          :: shortenedDt
   logical                          :: endRun

   !!!!!----- START INSERTION BY CODE GENERATOR
   integer              :: gpu_tf_hydro_nThreads
   integer              :: gpu_tf_hydro_nDistributorThreads
   integer              :: gpu_tf_hydro_nTilesPerPacket
   type(C_PTR)          :: gpu_tf_hydro_packet
   real(MILHOJA_REAL)   :: MH_dt
   integer(MILHOJA_INT) :: MH_ierr

   gpu_tf_hydro_packet = C_NULL_PTR 
   !!!!!----- END INSERTION BY CODE GENERATOR

   endRun = .FALSE.

   !!!!!----- START INSERTION BY CODE GENERATOR
   ! RPs are used directly by the Driver and therefore should be handled at
   ! this level rather than at the level of the code generated for use by
   ! the runtime.
   CALL RuntimeParameters_get("gpu_tf_hydro_nThreads", &
                               gpu_tf_hydro_nThreads)
   CALL RuntimeParameters_get("gpu_tf_hydro_nDistributorThreads", &
                               gpu_tf_hydro_nDistributorThreads)
   CALL RuntimeParameters_get("gpu_tf_hydro_nTilesPerPacket", &
                               gpu_tf_hydro_nTilesPerPacket)
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
      ! thread team configuration.  In this case, the GPU-only config.  The
      ! mapping should have been specified originally in the recipe.
      CALL Grid_fillGuardCells(CENTER, ALLDIR, &
                               doEos=.TRUE., &
                               maskSize=hy_gcMaskSize, &
                               mask=hy_gcMask, &
                               makeMaskConsistent=.TRUE., &
                               selectBlockType=LEAF, &
                               doLogMask=.TRUE.)

      MH_dt = REAL(dr_dt, kind=MILHOJA_REAL)
      MH_ierr = instantiate_hydro_advance_packet_C(MH_dt, &
                                                   gpu_tf_hydro_packet)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
      CALL Orchestration_executeTasks_Gpu(gpu_tf_hydro, &
                                          gpu_tf_hydro_nDistributorThreads, &
                                          gpu_tf_hydro_nThreads, &
                                          gpu_tf_hydro_nTilesPerPacket, &
                                          gpu_tf_hydro_packet)

      MH_ierr = delete_hydro_advance_packet_C(gpu_tf_hydro_packet)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
      gpu_tf_hydro_packet = C_NULL_PTR 
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

      ! The runtime backend's memory manager is presently very rudimentary.
      ! As a result, we need to ask it to reset itself at each step to avoid
      ! running out of memory on the GPU.
      CALL milhoja_runtime_reset(MH_ierr)
      CALL Orchestration_checkInternalError("Driver_evolveAll", MH_ierr)
   end do
   dr_nstep = min(dr_nstep, dr_nend)

   CALL Timers_stop("evolution")
   CALL Logfile_stamp('Exiting evolution loop', '[Driver_evolveAll]')

   call IO_outputFinal()

   CALL Timers_getSummary(MAX(0, dr_nstep - dr_nbegin + 1))
   CALL Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   CALL Logfile_close()
end subroutine Driver_evolveAll

