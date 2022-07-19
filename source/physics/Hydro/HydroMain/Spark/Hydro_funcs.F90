!!****i** source/physics/Hydro/HydroMain/Spark/Hydro_funcs
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
!!  Hydro_funcs
!!
!! DESCRIPTION
!!
!!  Holds functions frequently accessed by Hydro.F90 for Spark hydrodynamic solver.
!! 
!!***
!!Reorder directive used by FLASH with --index-reorder flag at setup
!!Reorder(4): hy_starState,Uin, U, hy_fl[xyz], hy_fluxBuf[XYZ]

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS


subroutine addFluxes(lev,blkLimits,weight,addFlux)
  !Store weighted fluxes, summed over RK stages, in temporary flux buffers.
  use Hydro_data, ONLY : hy_flx, hy_fly, hy_flz, hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ 
  
  implicit none
  
  integer,intent(IN)::lev
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
  real, intent(IN) :: weight
  logical, intent(IN) :: addFlux
  
  
  if (addFlux) then
     hy_fluxBufX = hy_fluxBufX+weight*hy_flx(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+0*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+0*K3D) 
     if (NDIM > 1) &   
          hy_fluxBufY = hy_fluxBufY+weight*hy_fly(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+0,&
 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1*K2D,&
 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+0*K3D)
     if (NDIM > 2) &
          hy_fluxBufZ = hy_fluxBufZ+weight*hy_flz(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+0,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+0*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1*K3D)
  else
     hy_fluxBufX = weight*hy_flx(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+0*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+0*K3D) 
     if (NDIM > 1) &   
          hy_fluxBufY = weight*hy_fly(:, blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+0,&
 blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1*K2D,&
 blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+0*K3D)
     if (NDIM > 2) &
          hy_fluxBufZ = weight*hy_flz(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+0,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+0*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1*K3D)
     
  endif
end subroutine addFluxes

!! Allocate variable size array holding local gravitational accelerations (depending on 
!! block size), fluxes, flux buffers, and save. 
!! If offloading to a device, send data and allocate on device only.
subroutine saveState(Uin,blkLimits,blklimitsGC)
  use Timers_interface, ONLY : Timers_start, Timers_stop
  use Hydro_data, ONLY : hya_starState, hya_tmpState
  implicit none
  
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  real, pointer :: Uin(:,:,:,:)
    real, pointer,dimension(:,:,:,:)::hy_starState,hy_tmpState
  integer ::  v,i1,i2,i3
  ! call Timers_start("Allocations")
  
  ! Allocate needed space on GPU if it is not already there
  
  hy_starState(1:NUNK_VARS,&
                blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_starState
  hy_tmpState(1:NUNK_VARS,&
                blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_tmpState
  ! update temp vars with solution data
  
  do i3=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do i2=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i1=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           do v=1,NUNK_VARS
              hy_starState(v,i1,i2,i3) = Uin(v,i1,i2,i3)
              hy_tmpState(v,i1,i2,i3) = Uin(v,i1,i2,i3)
           enddo
        end do
     end do
  end do
  nullify (hy_starState)
  nullify(hy_tmpState)
end subroutine saveState


subroutine updateState(Uin,blkLimits,blkLimitsGC)
  use Hydro_data, ONLY : hya_starState, hy_threadWithinBlock, hy_grav, hy_flx, hy_fly, hy_flz,&
                         hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ, hy_fluxCorrect, &
                         hy_smalldens, hy_smallE, hy_smallpres, &
                         hy_smallX, hy_cvisc, hy_del
  use Hydro_data, ONLY : hy_flat3d
  
  implicit none
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  real,dimension(:,:,:,:),pointer :: Uin
  real,dimension(:,:,:,:),pointer :: hy_starState

   hy_starState(1:NUNK_VARS,&
                blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
                blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
                blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS))=>hya_starState
#ifdef GPOT_VAR
  ! First reset GPOT_VAR.
  hy_starState(GPOT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       Uin(GPOT_VAR,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
#endif
  Uin(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)) = &
       hy_starState(:,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS),&
       blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS),&
       blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS))
  nullify(hy_starState)
end subroutine updateState

!! Set loop limits.  We include ngcell layers of guard zones
subroutine setLims(ngcell,blkLimits,limits)
  implicit none
  integer, intent(IN):: ngcell
  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: blkLimits
  integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: limits
  
  integer i
  
  
  limits = blkLimits
  do i=1,NDIM
     limits(LOW ,i) = blkLimits(LOW ,i) - ngcell
     limits(HIGH,i) = blkLimits(HIGH,i) + ngcell
  end do
end subroutine setLims

subroutine shockDetect(Uin,limits,blkLimitsGC)

  use Hydro_data,        only : hy_geometry, hy_tiny, hy_Vc

  implicit none


  !! ---- Argument List ----------------------------------
  real, dimension(:,:,:,:), pointer   :: Uin
  integer, intent(IN) :: limits(LOW:HIGH,MDIM)
  integer, intent(IN) :: blkLimitsGC(LOW:HIGH,MDIM)
  !! -----------------------------------------------------

  integer :: i,j,k
  logical :: SW1, SW2


  real :: divv,gradPx,gradPy,gradPz
  real :: minP,minC,beta,delta
  real :: localCfl,cflMax
  
  !necessary for argument for %getDataPtr()

#ifndef SHOK_VAR
  return
#endif

  ! Two parameters that can be adjusted to detect shocks
  ! with different strengths:
  ! (a) The lower the values the weaker shocks detected
  !     (lower the values to detect more shock regions)
  ! (b) The larger the values the stronger shocks detected
  !     (increase the values to detect less shock regions)
  beta = 0.5 !0.5 !10. ! gradP
  delta= 0.1  !0.1 !2. ! divV
!!$  beta  = 0.1 !0.1
!!$  delta = 0.01



  Uin(SHOK_VAR,:,:,:) = 0.
  

  !! Allocate a temporary cell-centered array for sound speed


  !! Compute sound speed
  do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
  do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
  do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
           hy_Vc(i,j,k) = sqrt(Uin(GAMC_VAR,i,j,k)*Uin(PRES_VAR,i,j,k)&
                /max(Uin(DENS_VAR,i,j,k),hy_tiny))
  end do
  end do
  end do

  do k=limits(LOW,KAXIS),limits(HIGH,KAXIS)
  do j=limits(LOW,JAXIS),limits(HIGH,JAXIS)
  do i=limits(LOW,IAXIS),limits(HIGH,IAXIS)

           ! initialize switch values
           SW1 = .false.
           SW2 = .false.

#if NDIM==1
           minP = minval(Uin(PRES_VAR,i-1:i+1,j,k))
           minC = minval(hy_Vc(i-1:i+1,j,k))
#endif
#if NDIM==2
           minP = minval(Uin(PRES_VAR,i-1:i+1,j-1:j+1,k))
           minC = minval(hy_Vc(i-1:i+1,j-1:j+1,k))
#endif
#if NDIM==3
           minP = minval(Uin(PRES_VAR,i-1:i+1,j-1:j+1,k-1:k+1))
           minC = minval(hy_Vc(i-1:i+1,j-1:j+1,k-1:k+1))
#endif
           !! We do not need to include non-Cartesian geom factors here.
           !! Undivided divV
           divv =        Uin(VELX_VAR,i+1,j,  k  ) - Uin(VELX_VAR,i-1,j,  k  )
#if NDIM > 1
           divv = divv + Uin(VELY_VAR,i,  j+1,k  ) - Uin(VELY_VAR,i,  j-1,k  )
#if NDIM == 3
           divv = divv + Uin(VELZ_VAR,i,  j,  k+1) - Uin(VELZ_VAR,i,  j,  k-1)
#endif
#endif
           divv = 0.5*divv

           !! Undivided grad pres
           gradPx = 0.5*(Uin(PRES_VAR,i+1,j,  k  ) - Uin(PRES_VAR,i-1,j,  k  ))
           gradPy = 0.
           gradPz = 0.
#if NDIM > 1
           gradPy = 0.5*(Uin(PRES_VAR,i,  j+1,k  ) - Uin(PRES_VAR,i,  j-1,k  ))
#if NDIM == 3
           gradPz = 0.5*(Uin(PRES_VAR,i,  j,  k+1) - Uin(PRES_VAR,i,  j,  k-1))
#endif
#endif
           if ( abs(gradPx)+abs(gradPy)+abs(gradPz) .ge. beta*minP ) then
              SW1 = .true.
           endif
           if (-delta*minC .ge. divv) then
              SW2 = .true.
           endif
           if (SW1 .and. SW2) then
              ! Set SHOCK_VAR to 1.0 if a shock is detected.
              ! One use is for a local hybrid method in the Hydro unit which
              ! applies (a diffusive) HLL solver when SHOK_VAR = 1.
              Uin(SHOK_VAR,i,j,k) = 1.
           endif !endif (SW1 .and. SW2) then
    end do
    end do
    end do
    

end subroutine shockDetect

!!Calculate divergence of the magnetic field.
subroutine calcDivB(Uin,hy_del,blkLimits)
  implicit none

  real, dimension(:,:,:,:), pointer   :: Uin
  real, dimension(MDIM),intent(IN) :: hy_del
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits


  real :: divB


#ifdef SPARK_GLM

  do k=blkLimits(LOW,KAXIS),blkLimits(HIGH,KAXIS)
  do j=blkLimits(LOW,JAXIS),blkLimits(HIGH,JAXIS)
  do i=blkLimits(LOW,IAXIS),blkLimits(HIGH,IAXIS)
           divB = 0.0
#if NDIM>1
           divB = (Uin(MAGX_VAR,i+1,j,k) - Uin(MAGX_VAR,i-1,j,k))&
                *0.5/hy_del(IAXIS)
           divB = divB + (Uin(MAGY_VAR,i,j+1,k) - Uin(MAGY_VAR,i,j-1,k))&
                *0.5/hy_del(JAXIS)
#if NDIM==3
           divB = divB + (Uin(MAGZ_VAR,i,j,k+1) - Uin(MAGZ_VAR,i,j,k-1))&
                *0.5/hy_del(KAXIS)
#endif
#endif
           Uin(DIVB_VAR,i,j,k) = divB
   end do
   end do
   end do
#endif
end subroutine calcDivB

