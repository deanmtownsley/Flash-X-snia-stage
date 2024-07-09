!!****if* source/Simulation/SimulationMain/IsentropicVortex/Simulation_init
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
!!  NAME
!!    Simulation_init
!!
!!  SYNOPSIS
!!    call Simulation_init( integer(in) : sim_meshMe )
!!
!!  DESCRIPTION
!!    Initialize all the data specified in Simulation_data
!!    Calls the RuntimeParameters_get to get the values of
!!    runtime parameters needed for the initial condition
!!    specification of IsentropicVortex problem
!!
!!    
!!  ARGUMENTS
!!    sim_meshMe    Current processor number
!!
!!
!!***

subroutine Simulation_init()

  use Simulation_data, ONLY : sim_gamma, sim_uAmbient, sim_vAmbient,&
       sim_vortexStrength, sim_xctrTrue, sim_yctrTrue, sim_nxSubint, &
       sim_nySubint,  sim_imax, sim_jmax, sim_imin, sim_jmin, &
       sim_imidDomain, sim_jmidDomain, sim_diDomain, sim_djDomain,&
       sim_tStarAmbient, sim_rbar, sim_constAmbient, sim_smlrho, sim_smallx,&
       sim_eosData, sim_eosMassFr, sim_meshMe
  use Driver_interface, ONLY : Driver_abort, Driver_getMype
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Eos_interface, ONLY : Eos


  implicit none

#include "constants.h"
#include "Simulation.h"
#include "Eos.h"

  real, save          :: rhoAmbient, pAmbient
  real          :: eAmbient, tAmbient,gm1i

  integer :: i, j, k, n
  integer :: ii, jj
  integer :: vecLen = 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real :: entropy, dst, dsd
  real :: pres, temp, dens, gamc, eint, entr, abar, zbar, ye
  call Driver_getMype(MESH_COMM,sim_meshMe)


  call RuntimeParameters_get( 'gamma', sim_gamma)
  call RuntimeParameters_get( 'u_ambient',sim_uAmbient  )
  call RuntimeParameters_get( 'v_ambient',sim_vAmbient  )
  call RuntimeParameters_get( 'vortex_strength', sim_vortexStrength  )
  call RuntimeParameters_get( 'xctr', sim_xctrTrue)
  call RuntimeParameters_get( 'yctr', sim_yctrTrue)
  call RuntimeParameters_get( 'nx_subint', sim_nxSubint)
  call RuntimeParameters_get( 'ny_subint', sim_nySubint)

  print *, sim_vortexStrength, sim_nxSubint, sim_nySubint

  call RuntimeParameters_get( 'smallx', sim_smallx)
  call RuntimeParameters_get( 'smlrho', sim_smlrho)
  
  call RuntimeParameters_get( 'xmax', sim_imax)
  call RuntimeParameters_get( 'xmin', sim_imin)
  call RuntimeParameters_get( 'ymax', sim_jmax)
  call RuntimeParameters_get( 'ymin', sim_jmin)

  call RuntimeParameters_get( 'rho_ambient', rhoAmbient  )
  call RuntimeParameters_get( 'p_ambient',   pAmbient  )
  
  call RuntimeParameters_get( 'diDomain', sim_diDomain)
  call RuntimeParameters_get( 'djDomain', sim_djDomain)
  call RuntimeParameters_get( 'imidDomain', sim_imidDomain)
  call RuntimeParameters_get( 'jmidDomain', sim_jmidDomain)
  
!   Some checking of the inputs:
  if ( NDIM == 1 ) then
     print *, 'Error: ndim = 1. Problem is designed for 2d or 3d.' 
     call Driver_abort('Error: ndim = 1; require ndim = 2 or 3.')
  endif
  
  if (rhoAmbient .LT. 1.e4*sim_smlrho) then
     print *, 'Error: ambient density is close to ', & 
          'cutoff density'
     print *, 'reset smlrho to be less than ', 1.e-4*rhoAmbient
     call Driver_abort('Error: ambient density too close to smlrho!')
  endif
!   Done checking.

!   Get t_ambient from rho_ambient and p_ambient.
  sim_eosMassFr=1.0
  dens=rhoAmbient
  pres=pAmbient
  call Eos(MODE_DENS_PRES,pres, temp, dens, gamc, eint, entr, abar, zbar, ye,sim_eosMassFr)
  tAmbient=temp
  eAmbient=eint
  
  
  !   Problem is actually specified in terms of a normalized temperature,
  !   which we will refer to as t_star. rbar is the scaling.
  sim_tStarAmbient = 1.0
  sim_rbar = sim_tStarAmbient/tAmbient
  gm1i = 1.0/(sim_gamma-1.0)
  sim_constAmbient = rhoAmbient/(tAmbient**gm1i)
  
  if (sim_meshMe .EQ. MASTER_PE) then
     
!!     call Logfile_stamp("initializing for Isentropic Vortex problem", 'run_init')
     print *, ' '
     print *, 'Flash-X: ', NDIM, 'dimensional vortex initialization'
     print *, ' '
     print *, 'Parameters read: '
     print *, ' '
     print *, 'gamma                 = ', sim_gamma
     print *, 'ambient density       = ', rhoAmbient
     print *, 'ambient pressure      = ', pAmbient
     print *, 'ambient x-velocity    = ', sim_uAmbient
     print *, 'ambient y-velocity    = ', sim_vAmbient
     print *, ' '
     print *, 'vortex_strength       = ', sim_vortexStrength
     print *, 'x center              = ', sim_xctrTrue
     print *, 'y center              = ', sim_yctrTrue
     print *, 'x subintervals        = ', sim_nxSubint
     print *, 'y subintervals        = ', sim_nySubint
     print *, ' '
     print *, 'Parameters computed : '
     print *, ' '
     print *, 'ambient temperature   = ', tAmbient
     print *, 'ambient int. energy   = ', eAmbient
     print *, 'gas constant          = ', sim_rbar
     print *, ' '
     print *, ' '
     
  endif
  
  
  return
end subroutine Simulation_init
