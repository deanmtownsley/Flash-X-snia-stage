!!****if* source/Simulation/SimulationMain/unitTest/Eos/Simulation_initBlock
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
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  call Simulation_initBlock(real,pointer :: solnData(:,:,:,:),
!!                            integer(IN)  :: tileDesc  )
!!
!!
!!
!! DESCRIPTION
!!   initializes a 3-d domain by setting with a density 
!!   gradient in x, temperature gradient in y, and composition gradient
!!   in z. An Eos call is made with mode = MODE_DENSITY_TEMP, and the
!!   resultant pressure and energy are saved. 
!!  
!!
!! 
!! ARGUMENTS
!!
!!  solnData  -        pointer to solution data
!!  tileDesc -        describes the block to initialize
!!
!!
!!***

!!REORDER(4): solnData
subroutine Simulation_initBlock(solnData,tileDesc)


  use Simulation_data, ONLY : sim_xmin,sim_xmax,sim_ymin,sim_ymax,&
                              sim_zmin,sim_zmax,sim_smallx,sim_smallE,&
                              sim_densMin, sim_tempMin, sim_yeMin, sim_presMin,&
                              sim_densMax, sim_tempMax, sim_yeMax, sim_presMax, &
                              sim_initialMass
  use Simulation_data, ONLY : sim_meshMe, sim_debug
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getGlobalIndexLimits
    
  use Grid_tile, ONLY : Grid_tile_t

  implicit none

#include "constants.h"
#include "Simulation.h"
  
  real,dimension(:,:,:,:),pointer :: solnData
  type(Grid_tile_t), intent(in) :: tileDesc

  integer :: iMax, jMax, kMax
  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC

  
  integer :: blockID
  integer :: i, j, k, n, even, put

  integer ::  npts_x, npts_y, npts_z
  real :: xpos, ypos, zpos
  real :: dx, dy, dz

  real ldensMax, ldensMin, dlogrho
  real ltempMax, ltempMin, dlogT
  real lpresMax, lpresMin, dlogP
  real yeMax, yeMin, dye
  real :: dens,temp,pres,ye

  
  real,allocatable,dimension(:) :: xCenter,yCenter,zCenter
  integer,dimension(MDIM) :: pos,globalInd
  integer,save :: callNo = 0
  

  call Grid_getGlobalIndexLimits(globalInd)

  blkLimits = tileDesc%limits
  blkLimitsGC = tileDesc%grownLimits

  if (sim_debug .AND. (sim_meshMe == MASTER_PE)) then
     callNo = callNo + 1
#ifdef FLASH_GRID_PARAMESH
     blockID = tileDesc%id
#else
     blockID = tileDesc % grid_index
#endif
     print*,'call Simulation_initBlock #',callNo,blockID
     print*,'globalInd  :',globalInd
     print*,'blkLimits  :',blkLimits
     print*,'blkLimitsGC:',blkLimitsGC
  end if

  allocate(xCenter( blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS) ))
  allocate(yCenter( blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS) ))
  allocate(zCenter( blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS) ))

  zCenter=0.0
  yCenter=0.0
#if NDIM == 3
  call Grid_getCellCoords(KAXIS, CENTER, tileDesc%level, &
                          blkLimitsGC(LOW,  :), &
                          blkLimitsGC(HIGH, :), &
                          zCenter)
#endif
#if NDIM >= 2
  call Grid_getCellCoords(JAXIS, CENTER, tileDesc%level, &
                          blkLimitsGC(LOW,  :), &
                          blkLimitsGC(HIGH, :), &
                          yCenter)
#endif
  call Grid_getCellCoords(IAXIS, CENTER, tileDesc%level, &
                          blkLimitsGC(LOW,  :), &
                          blkLimitsGC(HIGH, :), &
                          xCenter)

  npts_x = globalInd(IAXIS)
  npts_y = globalInd(JAXIS)
  npts_z = globalInd(KAXIS)
  
! compute the log interval for dens, temp, and xn
  dlogrho = (log10(sim_densMax) - log10(sim_densMin)) / float(npts_x)
  dlogT   = (log10(sim_tempMax) - log10(sim_tempMin)) / float(npts_y)
  dlogP   = (log10(sim_presMax) - log10(sim_presMin)) / float(npts_y)
  dye     = (sim_yeMax - sim_yeMin) / float(npts_z)

! compute the distance between coordinates -- ideally xmin=0, xmax=1
  dx = (sim_xmax - sim_xmin) / float(npts_x)
  dy = (sim_ymax - sim_ymin) / float(npts_y)
  dz = (sim_zmax - sim_zmin) / float(npts_z)

! since we will be making the density, etc. equally spaced in log, compute
! the log of the extrema
  ldensMax = log10(sim_densMax)
  ldensMin = log10(sim_densMin)
  
  ltempMax = log10(sim_tempMax)
  ltempMin = log10(sim_tempMin)
  
  lpresMax = log10(sim_presMax)
  lpresMin = log10(sim_presMin)
  
  yeMax    = sim_yeMax
  yeMin    = sim_yeMin

! now fill the master arrays
  do k = blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)

     zpos = (zCenter(k) - sim_zmin) / dz
     do j = blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
        
        ypos = (yCenter(j) - sim_ymin) / dy
        do i = blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           
           xpos = (xCenter(i) - sim_xmin) / dx
           ! compute the density, temperature, and composition
           dens = 10.e0**(ldensMin + xpos*dlogrho)
           temp = 10.e0**(ltempMin + ypos*dlogT)
           pres = 10.e0**(lpresMin + ypos*dlogP)
           ye   = yeMin + zpos*dye

           solnData(DENS_VAR,i,j,k)=dens
           solnData(CTMP_VAR,i,j,k)=temp
           solnData(TEMP_VAR,i,j,k)=temp
           solnData(CPRS_VAR,i,j,k)=pres
           solnData(YE_MSCALAR,i,j,k)=ye
           ! Put some small value in energies, so Eos_wrapped called at grid
           ! init will not have to floor then with smallE
           solnData(ENER_VAR,i,j,k)=sim_smallE*1.02
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k)=sim_smallE*1.02
#endif

           if (sim_debug) then
!!$              if (i==1 .AND. j==1 .AND. k==1) &
              print*,'i,j,k,dens,temp,pres:',i,j,k,dens,temp,pres
           end if
           
        enddo
     enddo
  enddo
  
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
  
  return
end subroutine Simulation_initBlock




