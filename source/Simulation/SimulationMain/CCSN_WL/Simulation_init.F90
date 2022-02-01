!!****if* source/Simulation/SimulationMain/CCSN_WL/Simulation_init
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
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!! DESCRIPTION
!!
!!  Initialization routine for CCSN problem.  Reads in 
!!  runtime parameters.  Reads 1D model file and maps
!!  variable names to FLASH unk names.
!!  Read-in data store in model_1d whose structure is defined
!!  by Simulation.h.
!!
!! ARGUMENTS
!!
!!   
!!
!! PARAMETERS
!!
!! NOTES
!!  
!!  This problem is described in, e.g.,
!!  Couch, S.M. 2013, ApJ, 765, 29
!!  Couch, S.M. 2013, ApJ, 775, 35
!!  Couch, S.M. & O'Connor, E.P. 2013, arXiv:1310.5728
!!
!!***
subroutine Simulation_init()
  use Simulation_data
  use Simulation_interface, ONLY : Simulation_mapIntToStr
  use Driver_interface, ONLY : Driver_getComm, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stampMessage
  use Grid_interface, ONLY : Grid_getGeometry

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Flashx_mpi.h"

  integer :: i

  integer, parameter :: max_stored_vars = 30
  real var_temp(max_stored_vars)

  character (len=256) :: current_line

  integer :: j, ipos, NUNK_VARS_stored, ierr
  integer :: var_key (NUNK_VARS)
  character (len=4) :: var_labels(max_stored_vars)

  call Driver_getComm(MESH_COMM, sim_meshComm)
  call Driver_getMype(MESH_COMM, sim_meshMe)

  call RuntimeParameters_get( 'model_file', model_file)
  call RuntimeParameters_get( 'smallx', sim_smallx)
  call RuntimeParameters_get( 'small', sim_small)
  call RuntimeParameters_get( 'vel_mult', sim_velMult)
  call RuntimeParameters_get( 'nsub', nsub)
  call RuntimeParameters_get( 'restart', sim_restart)

  do i = UNK_VARS_BEGIN, UNK_VARS_END
     call Simulation_mapIntToStr(i, unklabels(i), MAPBLOCK_UNK)
     call makeLowercase(unklabels(i))
  enddo

  ! open the file and read in the header 
  if ( (.not. sim_restart) .and. (sim_meshMe == MASTER_PE) ) then
    open(unit=2,file=model_file,status='old')
    print *, 'file opened'
    read (2,'(a80)') current_line
  
    ! read in the number of variables line
    read (2,'(a80)') current_line
    ipos = index(current_line,'=') + 1
    read (current_line(ipos:),*) nvar_stored
    print *,"read nvar_stored", nvar_stored
  
    if (NUNK_VARS .NE. nvar_stored) then
       print *, ' '
       print *, 'Warning: the number of variables stored in the'
       print *, 'input file is different than the number of'
       print *, 'variables in the current version of FLASH.'
       print *, ' '
       print *, 'The variables in the file that are also defined'
       print *, 'in FLASH will be read in.  Any missing variables'
       print *, 'will be initialized to zero'
       print *, ' '
    endif
  
    print *, "Vaiables in file:"

    do i = 1, nvar_stored
       read (2,'(a4)') var_labels(i)
       print *, var_labels(i)
       call makeLowercase(var_labels(i))
    enddo
  
    do j = 1, NUNK_VARS
       var_key(j) = NONEXISTENT
  
       do i = 1, nvar_stored
  
          if (unklabels(j) == var_labels(i)) then
             var_key(j) = i
          endif
  
       enddo
  
       if (var_key(j) == NONEXISTENT) then
          print *, 'Warning, variable: ', unklabels(j), & 
               ' not found in the input file.'
          print *, 'initializing ', unklabels(j), ' to 0'
          print *, ' '
       endif
  
    enddo
  
    xzn = 0.0
    model_1d = 0.0
    do i = 1, n1d_max
       read(2,*,end=11) xzn(i), (var_temp(j),j=1,nvar_stored)
  
       ! put these in order, so model1d_var always contains the same variables 
       ! in the same spots
       do j = 1, NUNK_VARS
          if (var_key(j) /= NONEXISTENT) then
             model_1d(i,j) = var_temp(var_key(j))
          else
             model_1d(i,j) = sim_small
          endif
  
       enddo
  
    enddo
  
11 close(unit=2)
  
    n1d_total = 0
  
    do while (xzn(n1d_total+1) .NE. 0.00)
       n1d_total = n1d_total + 1
    enddo
  
    print *, 'file read completed'
    print *, n1d_total, 'points read in'
    print *

  end if ! (.not. sim_restart)

  call MPI_Bcast (n1d_total,1,FLASH_INTEGER,MASTER_PE,sim_meshComm,ierr)
  call MPI_Bcast (xzn,n1d_max,FLASH_REAL,MASTER_PE,sim_meshComm,ierr)
  call MPI_Bcast (model_1d,n1d_max*NUNK_VARS,FLASH_REAL,MASTER_PE,sim_meshComm,ierr)

end subroutine Simulation_init
