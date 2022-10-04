!!****if* source/Simulation/SimulationMain/unitTest/Gravity/Poisson3/Driver_evolveAll
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
!!  Driver_evolveAll()
!!
!! DESCRIPTION
!!
!! This routine implements the Strang splitting scheme for time
!! advancement. A single step in the this driver 
!! includes two sweeps, the first one in order XYZ, and
!! the second one in order ZYX. This driver works with directionally
!! split operators only. The routine also controls the regridding of
!! the mesh if necessary and the simulation output.
!!
!!  
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***


#ifdef DEBUG_ALL
#define DEBUG_DRIVER
#endif
#define DEBUG_DRIVER

subroutine Driver_evolveAll()

  use Driver_data, ONLY: dr_globalMe, dr_nbegin, &
       dr_nend, dr_dt, dr_wallClockTimeLimit, &
       dr_simTime, dr_redshift, &
       dr_nstep, dr_dtOld, dr_dtNew, dr_restart, dr_elapsedWCTime
  use Driver_interface, ONLY : Driver_sourceTerms, Driver_computeDt, &
       Driver_getElapsedWCTime
  use Logfile_interface, ONLY : Logfile_stamp, Logfile_close
  use Timers_interface, ONLY : Timers_start, Timers_stop, &
    Timers_getSummary
  use Particles_interface, ONLY : Particles_advance, Particles_dump
  use Grid_interface, ONLY : Grid_getLocalNumBlks, &
    Grid_getListOfBlocks, Grid_updateRefinement
  use Hydro_interface, ONLY : Hydro
  use Gravity_interface, ONLY :  Gravity_potential, Gravity_unitTest
  !use IO_data, ONLY: io_justCheckpointed 
  use IO_interface, ONLY :IO_output,IO_outputFinal

  implicit none

#include "constants.h"
#include "Simulation.h"

  integer   :: localNumBlocks

  ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(3,2) :: strBuff
  character(len=15) :: numToStr

  ! for unit test
  character(len=4)            :: rank_str
  logical,          save      :: perfect = .true.
  character(len=20)           :: fileName
  integer,          parameter :: fileUnit = 2

  ! ------------ end of unitTest setup ---------------------------------------
  
  call Logfile_stamp( 'Entering evolution routine' , '[Driver_evolveAll]')

  

  call Timers_start("evolution")
  print*,' starting ',dr_nend, dr_nbegin
  if (dr_nend .GE. dr_nbegin) then


     if (dr_globalMe == MASTER_PE) then

        write (numToStr(1:), '(I10)') dr_nstep
        write (strBuff(1,1), "(A)") "n"
        write (strBuff(1,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_simTime
        write (strBuff(2,1), "(A)") "t"
        write (strBuff(2,2), "(A)") trim(adjustl(numToStr))
        
        write (numToStr(1:), "(1PE12.6)") dr_dt
        write (strBuff(3,1), "(A)") "dt"
        write (strBuff(3,2), "(A)") trim(adjustl(NumToStr))
        
        call Logfile_stamp( strBuff, 3, 2, "step")
     end if
     !--------------------------------------------------------------------
     !- Start Physics Sequence
     !----
  
     dr_simTime = dr_simTime + dr_dt

     call Timers_start("hydro")
#ifdef DEBUG_DRIVER
     print*,'going into hydro'
#endif

     call Hydro(dr_simTime, dr_dt, dr_dtOld, SWEEP_XYZ)

     call Timers_stop("hydro")

     
#ifdef DEBUG_DRIVER
     print*, 'return from Hydro/MHD timestep'
#endif

     call Timers_start("Particles_advance")
     call Particles_advance(dr_dtOld, dr_dt)
#ifdef DEBUG_DRIVER
     print*, 'return from Particles_advance '
#endif
     call Timers_stop("Particles_advance")     
     call Gravity_potential()
#ifdef DEBUG_DRIVER
     print*, 'return from Gravity_potential '
#endif

     dr_simTime = dr_simTime + dr_dt
     call Timers_start("hydro")
     call Hydro(dr_simTime, dr_dt, dr_dtOld, SWEEP_ZYX)
     call Timers_stop("hydro")



     call Timers_start("Particles_advance")
     call Particles_advance(dr_dt, dr_dt)
     call Timers_stop("Particles_advance")
     
     call Gravity_potential()

     !----
     !- End Physics Sequence
     !--------------------------------------------------------------------
  end if

  ! Gravity unitTest calculations-------------------------------------
  write(rank_str,"(I4.4)") dr_globalMe
  fileName = "unitTest_" // rank_str

  open(fileUnit,file=fileName)
  write(fileUnit,'("P",I0)') dr_globalMe

  call Timers_start("testing")
  call Gravity_unitTest(fileUnit,perfect)
  call Timers_stop("testing")

  call Timers_start("writing")
  if (perfect) then
    write(fileUnit,'("all results conformed with expected values.")')
  endif

  close(fileUnit)
  call Timers_stop("writing")

!! Eliminted all code beyond here, not needed for unit test -PMR
!! NO!  We'd actually like to SEE what was calculated. LBR
    !io_justCheckpointed = .false.
     
  ! ------------------------------- Gravity unitTest output

  call Timers_stop("evolution")

  call Logfile_stamp( 'Exiting evolution routine' , '[Driver_evolveAll]')

  call IO_outputFinal()

  call Timers_getSummary( dr_nstep)


  call Logfile_stamp( "FLASH run complete.", "LOGFILE_END")

  call Logfile_close()


  return
  
end subroutine Driver_evolveAll



