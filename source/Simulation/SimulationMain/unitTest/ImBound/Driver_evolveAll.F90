!!****if* source/Simulation/SimulationMain/unitTest/ImBound/Driver_evolveAll
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
!!
!!  Driver_evolveAll
!!
!! SYNOPSIS
!!
!!  call Driver_evolveAll()
!!
!! DESCRIPTION
!!
!!  Driver_evolveAll for incompFlow Simulations
!!
!!  DOC: Driver_evolveAll needs more explanation
!!
!! NOTES
!!
!!  variables that begin with "dr_" like, dr_globalMe or dr_dt, dr_beginStep
!!  are stored in the data fortran module for the Driver unit, Driver_data.
!!  The "dr_" is meant to indicate that the variable belongs to the Driver Unit.
!!  all other normally named variables i, j, etc are local variables.
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"
#include "FortranLangFeatures.fh"

#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif

subroutine Driver_evolveAll()

   use Driver_data, ONLY: dr_globalMe, dr_globalNumProcs, dr_nbegin, &
                          dr_nend, dr_dt, &
                          dr_tmax, dr_simTime, dr_redshift, &
                          dr_nstep, dr_dtOld, dr_dtNew, &
                          dr_simGeneration, &
                          dr_restart
   use Driver_interface, ONLY: Driver_sourceTerms, Driver_computeDt
   use Logfile_interface, ONLY: Logfile_stamp, Logfile_close
   use Timers_interface, ONLY: Timers_start, Timers_stop, &
                               Timers_getSummary
   use Particles_interface, ONLY: Particles_advance, Particles_dump

   use Grid_interface, ONLY: Grid_updateRefinement, Grid_setInterpValsGcell, &
                             Grid_fillGuardCells, Grid_getTileIterator, &
                             Grid_releaseTileIterator

   use Grid_iterator, ONLY: Grid_iterator_t

   use Grid_tile, ONLY: Grid_tile_t

   use IO_interface, ONLY: IO_output, IO_outputFinal

   use ImBound_interface, ONLY: ImBound_mapToGrid, ImBound_redistanceLS, &
                                ImBound_advance, ImBound_reInitGridVars, &
                                ImBound_getBodyPtr, ImBound_releaseBodyPtr, &
                                ImBound_skipBox

   use SolidMechanics_interface, ONLY: SolidMechanics_updateBodyForce

   use ImBound_data, ONLY: ib_lsIt, ib_numBodies, ib_enableSelectiveMapping

   use ImBound_type, ONLY: ImBound_type_t

   implicit none

   ! for logfile output
   character(len=MAX_STRING_LENGTH), dimension(3, 2) :: strBuff
   character(len=15) :: numToStr

   logical :: gridChanged
   logical :: endRunPl !Should we end our run on this iteration, based on conditions detected by the IO unit?
   logical :: endRun !Should we end our run on this iteration, based on conditions detected by the IO unit?
   logical :: endRunWallClock !Should we end our run on this iteration, based on wall clock time?

   !---------Variables for unitTest-----------
   character(len=20) :: fileName
   integer, parameter        :: fileUnit = 2
   integer, dimension(4) :: prNum
   integer :: temp, i
   logical :: gcMask(NUNK_VARS + NDIM*NFACE_VARS)
   type(Grid_iterator_t) :: itor
   type(Grid_tile_t) :: tileDesc
   type(ImBound_type_t), pointer :: bodyInfo
   integer :: iteration, ibd
   logical :: skipBox = .FALSE.

   nullify (bodyInfo)

   endRunPl = .false.
   endRun = .false.

   call Logfile_stamp('Entering evolution loop', '[Driver_evolveAll]')
   call Timers_start("evolution")

   ! Initial Timestep:
   ! backup needed old
   dr_dtOld = dr_dt

   ! calculate new
   call Driver_computeDt(dr_nbegin, dr_nstep, &
                         dr_simTime, dr_dtOld, dr_dtNew)
   ! store new
   dr_dt = dr_dtNew

   !------------------------------------------------------------
   do ibd = 1, ib_numBodies
      !------------------------------------------------------------
      call ImBound_getBodyPtr(bodyInfo, ibd)
      call Grid_getTileIterator(itor, nodetype=LEAF)
      !------------------------------------------------------------
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call ImBound_mapToGrid(tileDesc, bodyInfo)
         call itor%next()
      end do
      !------------------------------------------------------------
      call Grid_releaseTileIterator(itor)
      call ImBound_releaseBodyPtr(bodyInfo, ibd)
      !------------------------------------------------------------
   end do

   gcMask = .FALSE.
   gcMask(LMDA_VAR) = .TRUE.
   call Grid_fillGuardCells(CENTER, ALLDIR, &
                            maskSize=NUNK_VARS, mask=gcMask)

   do dr_nstep = dr_nBegin, dr_nend

      if (dr_globalMe == MASTER_PE) then

         write (numToStr(1:), '(I10)') dr_nstep
         write (strBuff(1, 1), "(A)") "n"
         write (strBuff(1, 2), "(A)") trim(adjustl(numToStr))

         write (numToStr(1:), "(1PE12.6)") dr_simTime
         write (strBuff(2, 1), "(A)") "t"
         write (strBuff(2, 2), "(A)") trim(adjustl(numToStr))

         write (numToStr(1:), "(1PE12.6)") dr_dt
         write (strBuff(3, 1), "(A)") "dt"
         write (strBuff(3, 2), "(A)") trim(adjustl(NumToStr))

         call Logfile_stamp(strBuff, 3, 2, "step")
      end if

      !------------------------------------------------------------
      !- Start Physics Sequence
      !------------------------------------------------------------
      dr_simTime = dr_simTime + dr_dt
      dr_simGeneration = 0
      !------------------------------------------------------------

      !------------------------------------------------------------
      call Grid_getTileIterator(itor, nodetype=LEAF)
      do while (itor%isValid())
         call itor%currentTile(tileDesc)
         !---------------------------------------------------------
         call ImBound_reInitGridVars(tileDesc)
         !---------------------------------------------------------
         call itor%next()
      end do
      call Grid_releaseTileIterator(itor)
      !------------------------------------------------------------

      !------------------------------------------------------------
      do ibd = 1, ib_numBodies
         !------------------------------------------------------------
         call ImBound_getBodyPtr(bodyInfo, ibd)
         call Grid_getTileIterator(itor, nodetype=LEAF)
         !------------------------------------------------------------
         do while (itor%isValid())
            call itor%currentTile(tileDesc)
            !---------------------------------------------------------
            call ImBound_skipBox(tileDesc, bodyInfo, skipBox)
            !---------------------------------------------------------
            if (skipBox) then
               !------------------------------------------------------
               call itor%next()
               cycle
            else
               !------------------------------------------------------
               call ImBound_mapToGrid(tileDesc, bodyInfo)
               call itor%next()
            end if
         end do
         !------------------------------------------------------------
         call Grid_releaseTileIterator(itor)
         call ImBound_releaseBodyPtr(bodyInfo, ibd)
         !------------------------------------------------------------
      end do

      gcMask = .FALSE.
      gcMask(LMDA_VAR) = .TRUE.
      call Grid_fillGuardCells(CENTER, ALLDIR, &
                               maskSize=NUNK_VARS, mask=gcMask)

      ! Apply redistancing procedure
      !------------------------------------------------------------
      do iteration = 1, ib_lsIt

         call Grid_getTileIterator(itor, nodetype=LEAF)
         do while (itor%isValid())
            call itor%currentTile(tileDesc)
            !------------------------------------------------------
            call ImBound_redistanceLS(tileDesc, iteration)
            !------------------------------------------------------
            call itor%next()
         end do
         call Grid_releaseTileIterator(itor)

         ! Fill GuardCells for level set function
         gcMask = .FALSE.
         gcMask(LMDA_VAR) = .TRUE.
         call Grid_fillGuardCells(CENTER, ALLDIR, &
                                  maskSize=NUNK_VARS, mask=gcMask)
      end do

      do ibd = 1, ib_numBodies
         call ImBound_getBodyPtr(bodyInfo, ibd)
         call SolidMechanics_updateBodyForce(bodyInfo, dr_simTime, dr_dt)
         call ImBound_advance(bodyInfo, dr_simTime, dr_dt)
         call ImBound_releaseBodyPtr(bodyInfo, ibd)
      end do

      !output a plotfile before the grid changes
      call Timers_start("IO_output")

      call IO_output(dr_simTime, &
                     dr_dt, dr_nstep + 1, dr_nbegin, endRunPl, PLOTFILE_AND_PARTICLEFILE)
      call Timers_stop("IO_output")

      ! Update grid and notify changes to other units
      call Grid_updateRefinement(dr_nstep, dr_simTime, gridChanged)

      if (dr_globalMe .eq. MASTER_PE) then
         write (*, *) ' '
         write (*, '(I6,A,g16.8,A,g16.8)') dr_nstep, &
            ', TimeStep= ', dr_dt, ', SimTime= ', dr_simTime
      end if

      if (dr_globalMe .eq. MASTER_PE) &
         write (*, *) '###############################################################################'

      ! Compute next step dt:
      ! backup needed old
      dr_dtOld = dr_dt

      ! calculate new
      call Driver_computeDt(dr_nbegin, dr_nstep, &
                            dr_simTime, dr_dtOld, dr_dtNew)
      ! store new
      dr_dt = dr_dtNew

      call Timers_start("io")
      call IO_output(dr_simTime, dr_dt, dr_nstep + 1, dr_nbegin, endRun, &
                     CHECKPOINT_FILE_ONLY)
      call Timers_stop("io")
      endRun = (endRunPl .OR. endRun)

     !!*****************************************************************************
     !!  Evolution Loop -- check termination conditions
     !!*****************************************************************************

      !Exit if a .dump_restart or .kill was found during the last step
      if (endRun) exit

     !! the simulation ends before nend iterations if
     !!  (i)   the simulation time is greater than the maximum time (tmax)
     !!  (ii)  the redshift falls below the minimum redshift
     !!        (also called zfinal)
     !!  (iii) the wall clock time is greater than the maximum
     !!        (wall_clock_time_max)

      if (dr_simTime >= dr_tmax) then
         if (dr_globalMe == MASTER_PE) then
            print *, "exiting: reached max SimTime"
         end if
         exit
      end if

      call dr_wallClockLimitExceeded(endRunWallClock)
      if (endRunWallClock) then
         if (dr_globalMe == MASTER_PE) then
            print *, "exiting: reached max wall clock time"
         end if
         exit
      end if

   end do
   dr_nstep = min(dr_nstep, dr_nend)

  !!******************************************************************************
  !! End of Evolution Loop
  !!******************************************************************************

   call Timers_stop("evolution")
   call Logfile_stamp('Exiting evolution loop', '[Driver_evolveAll]')
   if (.NOT. endRun) call IO_outputFinal()
   call Timers_getSummary(max(0, dr_nstep - dr_nbegin + 1))
   call Logfile_stamp("FLASH run complete.", "LOGFILE_END")
   call Logfile_close()

   temp = dr_globalMe
   do i = 1, 4
      prNum(i) = mod(temp, 10)
      temp = temp/10
   end do

   !-------------------------------------------------------------------------------

   return

end subroutine Driver_evolveAll
