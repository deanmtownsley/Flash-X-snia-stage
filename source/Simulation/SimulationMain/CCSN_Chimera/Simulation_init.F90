!!****if* source/Simulation/SimulationMain/CCSN_Chimera/Simulation_init
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  Simulation_init
!!
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
!!  Sandoval et al. 2021, ApJ, 921 113
!!
!!***
subroutine Simulation_init()
  use Simulation_data
  use Driver_interface, ONLY : Driver_abort, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_interface, ONLY : Grid_getGeometry, Grid_getDomainBoundBox
  !use gr_mpoleData, ONLY : gr_point_mass=>point_mass

  use chimera_model_module
  use model_interp_module

  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Flashx_mpi.h"

  integer :: i, j, k, l
  integer :: meshGeom

  real, dimension(LOW:HIGH,MDIM) :: boundBox
  real :: mass_chim, vol_chim, mass_prog, vol_prog
  real :: r, theta, phi, x, y, z, rcyl, zcyl, gr_max_r, point_mass
  real :: rlo, rhi, dr, dvol, dvolr, domega, domega_exclude, tmp1, tmp2, tmp3
  integer :: irho_inner

  call Driver_getMype(MESH_COMM, sim_meshMe)

  call Grid_getGeometry(meshGeom)
  call Grid_getDomainBoundBox(boundBox)

  call RuntimeParameters_get( 'chimera_model_file', chimera_model_file)
  call RuntimeParameters_get( 'progenitor_model_file', progenitor_model_file)
  call RuntimeParameters_get( 'smlrho',   sim_smlrho)
  call RuntimeParameters_get( 'smallt',   sim_smallt)
  call RuntimeParameters_get( 'smallx',   sim_smallx)
  call RuntimeParameters_get( 'restart', sim_restart)
  call RuntimeParameters_get( 'max_r', sim_max_r)
  call RuntimeParameters_get( 'r_inner', sim_r_inner)
  call RuntimeParameters_get( 'bdry_dens', sim_bdry_dens)
  call RuntimeParameters_get( 'bdry_pres', sim_bdry_pres)
  call RuntimeParameters_get( 'bdry_temp', sim_bdry_temp)

  ! load chimera model
  if ( len_trim(chimera_model_file) > 0 ) then
     call open_chimera_file(chimera_model_file)
     call read_chimera_file
  end if

  ! load progenitor model
  if ( len_trim(progenitor_model_file) > 0 ) then
    call read_progenitor_file(progenitor_model_file)
  end if

  point_mass = 0.0
  do i = 1, imax_chim+1
     if ( x_e_chim(i) < sim_r_inner ) then
        if ( x_e_chim(i+1) <= sim_r_inner ) then
           point_mass = point_mass + sum( dmass_e_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) )
        else
           dr = x_e_chim(i+1) - sim_r_inner
           dvolr = dr * ( sim_r_inner * x_e_chim(i+1) + dr * dr * (1.0/3.0) )
           point_mass = point_mass + ( dvolx_e_chim(i) - dvolr ) &
           &                         * sum( rho_c_chim(i,jmin_chim:jmax_chim,kmin_chim:kmax_chim) &
           &                                * domega_chim(jmin_chim:jmax_chim,kmin_chim:kmax_chim) )
        end if
     end if
  end do

  if ( sim_meshMe == MASTER_PE ) then
     write(*,*) 'point_mass=',point_mass
     write(*,*) 'adding this value to user-supplied value'
  end if

  !Excised CHIMERA point mass will be added to the point mass runtime parameter used in gravity routines
  !gr_point_mass = gr_point_mass + sim_pointMass
  sim_pointMass = point_mass

  if (sim_restart) then
     call Simulation_finalize()
  endif     

  contains

!!! THIS IS OLD READ IN OF KEPLER FILE
!!!----------------------------------------------------------------------------
     subroutine read_progenitor_file(progenitor_model_file)

        use Simulation_interface, ONLY : Simulation_mapIntToStr

        implicit none

        ! input variables
        character (*), intent(in) :: progenitor_model_file

        ! local variables
        integer, parameter :: max_stored_vars = 180
        real :: var_temp(max_stored_vars)
        integer :: ipos, NUNK_VARS_stored
        integer :: var_key(NUNK_VARS)
        character (len=4) :: var_labels(max_stored_vars)
        character (len=256) :: current_line
        integer :: prog_unit
        real :: dr
        integer :: i, j, k, ierr

        do i = UNK_VARS_BEGIN, UNK_VARS_END
           call Simulation_mapIntToStr(i, unklabels(i), MAPBLOCK_UNK)
           call makeLowercase(unklabels(i))
        enddo

        ! open the file and read in the header 
        open (newunit=prog_unit,file=trim(progenitor_model_file),status='old')
        read (prog_unit,'(a256)') current_line

        if (sim_meshMe == MASTER_PE) then
           print *, 'file opened'
        end if

        ! read in the number of variables line
        read (prog_unit,'(a256)') current_line
        ipos = index(current_line,'=') + 1
        read (current_line(ipos:),*) nvar_stored
        if (sim_meshMe == MASTER_PE) then
           print *,"read nvar_stored", nvar_stored
        end if

        if (NUNK_VARS /= nvar_stored .AND. sim_meshMe == MASTER_PE) then
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

        if (sim_meshMe == MASTER_PE) then
           print *, "Vaiables in file:"
        endif

        do i = 1, nvar_stored
           read (prog_unit,'(a4)') var_labels(i)
           if (sim_meshMe == MASTER_PE) then
              print *, var_labels(i)
           end if
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
              if(sim_meshMe == MASTER_PE) then
                 print *, 'Warning, variable: ', unklabels(j), ' not found in the input file.'
                 print *, 'initializing ', unklabels(j), ' to 0'
                 print *, ' '
              endif
           endif
        enddo

        do i = 1, n1d_max
           read(prog_unit,*,iostat=ierr) xzn(i), (var_temp(j),j=1,nvar_stored)
           if ( ierr /= 0 ) exit
           ! put these in order, so model1d_var always contains the same variables in the same spots
           do j = 1, NUNK_VARS
              if (var_key(j) /= NONEXISTENT) then
                 model_1d(i,j) = var_temp(var_key(j))
              else
                 model_1d(i,j) = 0.0
              endif
           enddo
        enddo
        close(prog_unit)

        n1d_total = 0
        do while (xzn(n1d_total+1) /= 0.0)
           n1d_total = n1d_total + 1
        enddo

        if (sim_meshMe == MASTER_PE) then
           print *, 'file read completed'
           print *, n1d_total, 'points read in'
        endif

        volxzn(1) = (1.0/3.0) * xzn(1)**3
        do j = 1, n1d_total
           dr          = xzn(j+1) - xzn(j)
           dvolxzn(j)  = dr * ( xzn(j) * xzn(j+1) + dr * dr * (1.0/3.0) )
           volxzn(j+1) = volxzn(j) + dvolxzn(j)
        end do

        xzn_ctr(1) = xzn(1)/2.0
        do k = 2, n1d_total
           xzn_ctr(k)  = ( xzn(k-1) + xzn(k) ) / 2.0
        enddo

  end subroutine read_progenitor_file
!!!----------------------------------------------------------------------------

end subroutine Simulation_init
