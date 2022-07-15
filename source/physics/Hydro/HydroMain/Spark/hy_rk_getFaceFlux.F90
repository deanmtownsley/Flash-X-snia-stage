!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_getFaceFlux
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
!!
!!  hy_rk_getFaceFlux
!!
!!  SYNOPSIS
!!
!!  call hy_rk_getFaceFlux ( integer(IN) :: blklimits(:,:),
!!                           integer(IN) :: blklimitsGC(:,:),
!!                           real(IN)    :: hy_del(:),
!!                           integer(IN) :: limits(:,:))
!!
!!  DESCRIPTION
!!  Initially stores grid data by 'pencils' (see subroutine setPencil() below),
!!  calls outside hy_reconstruction subroutine, calls outisde Riemann solver, and 
!!  ultimately yields face fluxes for relevant directions.  These fluxes are lastly 
!!  saved to the grid multifabs (AMReX) for later access.
!!
!!  ARGUMENTS
!!     blklimits,blkLimitsGC are the bounds of the block/tile
!!     hy_del are dx,dy,dz
!!     limits   sets limits of the fluxes.
!!
!!     Note this is a modified form typically yielded
!!     by blockDesc%limits b/c stage 1 of RK integration accounts for NSTENCIL 
!!     number of guard cells beyond the block interior.  Stage 2 only updates the 
!!    interior cells.
!!
!!***

!!Reorder(4):hy_starState,hy_fl[xyz] 
subroutine hy_rk_getFaceFlux (blklimits,blkLimitsGC, limits)

  use Hydro_data, ONLY : hy_threadWithinBlock, &
       hy_starState, hy_grav, hy_flattening, hy_flx, hy_fly, hy_flz
  use Hydro_data, ONLY : hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_cvisc
  use Hydro_data, ONLY : hy_rope, hy_uPlus, hy_uMinus, hy_flat, hy_grv, hy_shck, hy_flux, &
  hy_flat3d
  use Hydro_data, ONLY : hy_del, hy_dlim
  use hy_rk_interface, ONLY : hy_reconstruct
  use Timers_interface, ONLY : Timers_start, Timers_stop

  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS


  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits, blkLimits, blkLimitsGC
  real, pointer,dimension(:,:,:,:) :: leftState, rightState
  integer :: i1,i2,i3, n, g, v, dir, ierr,i,j,k,i_s,j_s,k_s
  integer, dimension(3) :: gCells
  integer, dimension(LOW:HIGH,MDIM) :: klim,dirbnds
  character(len = 2) :: dir_str
  real :: cvisc, VenerLo, VenerHi, accelM, accelP,dx
  integer :: s	 
  real :: spcSumInv
  real, pointer :: spc(:)
  
  logical :: inShock
  real, dimension(HY_NUM_VARS) :: VL, VR
  real :: rope(5)

  !$omp target data map(to: dir, klim,hy_dlim,gCells)
  if (hy_flattening) then
     call flattening(limits)
  else
     ! call Timers_start("hy_flat3d")
     !$omp target teams distribute parallel do collapse(3) &
     !$omp shared(blkLimitsGC,hy_flat3d) private(i1,i2,i3) default(none) map(to:blkLimitsGC)
     !! TODO: Set this once for both rk steps.
     do k=blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j=blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
     do i=blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                  hy_flat3d(i,j,k) = 1.0
     end do
     end do
     end do
     ! call Timers_stop("hy_flat3d")
  end if

  
   !  Begin loop over zones
  do dir = 1, NDIM
     select case(dir)
     case (IAXIS)
     hy_dlim(:,:) = limits(:,:)
     gCells(1) = (limits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
     gCells(2) = (limits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
     gCells(3) = (limits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))
     dir_str = "_x"
     case (JAXIS)
     hy_dlim(:,1) = limits(:,JAXIS)
     hy_dlim(:,2) = limits(:,IAXIS)
     hy_dlim(:,3) = limits(:,KAXIS)
     gCells(1) = (limits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
     gCells(2) = (limits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
     gCells(3) = (limits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))
     dir_str = "_y"
     case (KAXIS)
     hy_dlim(:,1) = limits(:,KAXIS)
     hy_dlim(:,2) = limits(:,IAXIS)
     hy_dlim(:,3) = limits(:,JAXIS)
     gCells(1) = (limits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))
     gCells(2) = (limits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
     gCells(3) = (limits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
     dir_str = "_z"
     end select
     

     
    
           !Define appropriate changing indices
           klim(LOW,:)=gCells(:)+1
           klim(LOW,1)=gCells(1)
           klim(HIGH,:)=1 + gCells(:) + hy_dlim(HIGH,:) - hy_dlim(LOW,:)
           klim(HIGH,1)=klim(HIGH,1)+1
           !$omp target update to(klim, dir, hy_dlim,gCells)
           !$omp target teams distribute parallel do collapse(3) default(none) &
           !$omp shared(dir,hy_starState,hy_rope,hy_grv,hy_grav,hy_shck,hy_flat,&
           !$omp hy_flat3d,klim,gCells,hy_dlim) private(i,j,k,n,i_s, j_s, k_s)
           do k = hy_dlim(LOW,3) - gCells(3), hy_dlim(HIGH,3) + gCells(3)
           do j = hy_dlim(LOW,2) - gCells(2), hy_dlim(HIGH,2) + gCells(2)
           do i = hy_dlim(LOW,1) - gCells(1), hy_dlim(HIGH,1) + gCells(1)
           i_s = 1 + i - hy_dlim(LOW,1) + gCells(1)
           j_s = 1 + j - hy_dlim(LOW,2) + gCells(2)
           k_s = 1 + k - hy_dlim(LOW,3) + gCells(3)
           if (dir == IAXIS) then
              hy_rope(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,i,j,k)
              hy_rope(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,i,j,k)
              hy_rope(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,i,j,k)
              hy_rope(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,i,j,k)
              hy_rope(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,i,j,k)
              hy_rope(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,i,j,k)
              hy_rope(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,i,j,k)*hy_starState(EINT_VAR,i,j,k)
#ifdef SPARK_GLM
              hy_rope(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,i,j,k)
              hy_rope(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,i,j,k)
              hy_rope(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,i,j,k)
              hy_rope(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,i,j,k)
#endif
#if NSPECIES+NMASS_SCALARS>0
              do n=SPECIES_BEGIN, MASS_SCALARS_END
              hy_rope(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,i,j,k)
              enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
              hy_grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,i,j,k)
#else
              hy_grv(i_s,j_s,k_s) = hy_grav(IAXIS,i,j,k)
#endif
#endif
#ifdef SHOK_VAR
              hy_shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,i,j,k)
#else
              hy_shck(i_s,j_s,k_s) = 0.0
#endif
              hy_flat(i_s,j_s,k_s) = hy_flat3d(i,j,k)
           else if (dir == JAXIS) then
              hy_rope(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,i,k)
              hy_rope(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,j,i,k)
              hy_rope(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,j,i,k)
              hy_rope(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,j,i,k)
              hy_rope(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,j,i,k)
              hy_rope(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,j,i,k)
              hy_rope(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,i,k)*hy_starState(EINT_VAR,j,i,k)
#ifdef SPARK_GLM
              hy_rope(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,j,i,k)
              hy_rope(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,j,i,k)
              hy_rope(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,j,i,k)
              hy_rope(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,j,i,k)
#endif
#if NSPECIES+NMASS_SCALARS>0
              do n=SPECIES_BEGIN, MASS_SCALARS_END
              hy_rope(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,j,i,k)
              enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
              hy_grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,j,i,k)
#else
              hy_grv(i_s,j_s,k_s) = hy_grav(IAXIS,j,i,k)
#endif
#endif
#ifdef SHOK_VAR
              hy_shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,j,i,k)
#else
              hy_shck(i_s,j_s,k_s) = 0.0
#endif
              hy_flat(i_s,j_s,k_s) = hy_flat3d(j,i,k)
           else if (dir == KAXIS) then
              hy_rope(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,k,i)
              hy_rope(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,j,k,i)
              hy_rope(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,j,k,i)
              hy_rope(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,j,k,i)
              hy_rope(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,j,k,i)
              hy_rope(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,j,k,i)
              hy_rope(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,k,i)*hy_starState(EINT_VAR,j,k,i)
#ifdef SPARK_GLM
              hy_rope(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,j,k,i)
              hy_rope(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,j,k,i)
              hy_rope(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,j,k,i)
              hy_rope(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,j,k,i)
#endif
#if NSPECIES+NMASS_SCALARS>0
              do n=SPECIES_BEGIN, MASS_SCALARS_END
              hy_rope(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,j,k,i)
              enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
              hy_grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,j,k,i)
#else
              hy_grv(i_s,j_s,k_s) = hy_grav(IAXIS,j,k,i)
#endif
#endif
#ifdef SHOK_VAR
              hy_shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,j,k,i)
#else
              hy_shck(i_s,j_s,k_s) = 0.0
#endif
              hy_flat(i_s,j_s,k_s) = hy_flat3d(j,k,i)
           endif
           end do
           end do
           end do  

           !$omp target teams distribute parallel do collapse(2) & ! This collapse 2 is because there is a data dependency
           !$omp private(i1,i2,i3) shared(hy_dlim,klim)

           ! call Timers_start("recon"//dir_str)
           
           !$omp target teams distribute parallel do collapse(4) &
           !$omp  private(i1,i2,i3,v) shared(klim,hy_uPlus, hy_uMinus, hy_rope, hy_flat) default(none)
           do i3=klim(LOW,KAXIS),klim(HIGH,KAXIS)
           do i2=klim(LOW,JAXIS),klim(HIGH,JAXIS)
           do i1=klim(LOW,IAXIS),klim(HIGH,IAXIS)      
              do v=1,NRECON
               call hy_reconstruct(v, hy_rope(v,i1-2:i1+2,i2,i3),hy_uPlus(v,i1,i2,i3),&
                    & hy_uMinus(v,i1,i2,i3), hy_flat(i1,i2,i3))
               enddo
           end do
           end do
           end do
               

           !$omp target teams distribute parallel do collapse(3) &
           !$omp private(i1,i2,i3,v) shared(dir,klim)
           do i3=klim(LOW,KAXIS),klim(HIGH,KAXIS)
           do i2=klim(LOW,JAXIS),klim(HIGH,JAXIS)
           do i1=klim(LOW,IAXIS),klim(HIGH,IAXIS)

           leftState => hy_uPlus
           rightState => hy_uMinus
           leftState(HY_DENS ,i1,i2,i3) = max(hy_smalldens, leftState(HY_DENS ,i1,i2,i3))
           leftState(HY_PRES ,i1,i2,i3) = max(hy_smallpres, leftState(HY_PRES ,i1,i2,i3))
#if NSPECIES>0
           ! Limit and renormalize the species.
           spc => leftState(HY_NUM_VARS+1:HY_NUM_VARS+NSPECIES ,i1,i2,i3)
           do s = 1, NSPECIES
           spc(s) = max(hy_smallX,min(1.0,spc(s)))
           end do
           spcSumInv = 1./sum(spc(1:NSPECIES))
           spc = spc*spcSumInv
#endif
           rightState(HY_DENS ,i1,i2,i3) = max(hy_smalldens, rightState(HY_DENS ,i1,i2,i3))
           rightState(HY_PRES ,i1,i2,i3) = max(hy_smallpres, rightState(HY_PRES ,i1,i2,i3))
#if NSPECIES>0
           ! Limit and renormalize the species.
           spc => rightState(HY_NUM_VARS+1:HY_NUM_VARS+NSPECIES ,i1,i2,i3)
           do s = 1, NSPECIES
           spc(s) = max(hy_smallX,min(1.0,spc(s)))
           end do
           spcSumInv = 1./sum(spc(1:NSPECIES))
           spc = spc*spcSumInv
#endif
           nullify(leftState)
           nullify(rightState)

           end do
           end do
           end do 

           
           !$omp target teams distribute parallel do collapse(2) & ! This collapse 2 is because there is a data dependency
           !$omp private(i1,i2,i3) shared(dir,klim)
           
              
           ! Check for shocks in the zones involved in flux calculation
           

           ! Now call the Riemann solver to compute fluxes
           !!klim(LOW,1)=klim(LOW,1)+1
           !$omp target teams distribute parallel do collapse(3) &
           !$omp private(i1,i2,i3) shared(dir,klim,hy_flux) default(none)
           do i3=klim(LOW,KAXIS),klim(HIGH,KAXIS)
           do i2=klim(LOW,JAXIS),klim(HIGH,JAXIS)
              do i1=klim(LOW,IAXIS),klim(HIGH,IAXIS)
                   VL = hy_uPlus(1:HY_NUM_VARS,i1-1,i2,i3)
                   VR = hy_uMinus(1:HY_NUM_VARS,i1,i2,i3)
                   inShock = any(hy_shck(i1-1:i1,i2,i3) /= 0.0)
                   call hy_riemann(dir,VL,VR,inShock,hy_flux(1:HY_NUM_FLUX,i1,i2,i3))
                end do
           end do
           end do
           
           !$omp target teams distribute parallel do collapse(3) &
           !$omp private(i1,i2,i3,cvisc,venerLo,venerHies) shared(hy_cvisc,hy_rope,dir,klim,hy_flux) default(none)
           do i3=klim(LOW,KAXIS),klim(HIGH,KAXIS)
           do i2=klim(LOW,JAXIS),klim(HIGH,JAXIS)
           do i1=klim(LOW,IAXIS),klim(HIGH,IAXIS)      
           ! Add artificial viscosity for strong-shock capturing
           cvisc = hy_cvisc*max(-(hy_rope(HY_VELX+dir-1,i1,i2,i3) - hy_rope(HY_VELX+dir-1,i1-1,i2,i3)),0.)
           
           ! Construct minus and plus TOTAL energy densities
           VenerLo = hy_rope(HY_DENS,i1-1,i2,i3)*0.5*(dot_product(hy_rope(HY_VELX:HY_VELZ,i1-1,i2,i3),hy_rope(HY_VELX:HY_VELZ,i1-1,i2,i3)))&
           + hy_rope(HY_RHOE,i1-1,i2,i3)
           VenerHi = hy_rope(HY_DENS,i1,i2,i3)*0.5*(dot_product(hy_rope(HY_VELX:HY_VELZ,i1,i2,i3),hy_rope(HY_VELX:HY_VELZ,i1,i2,i3)))&
           + hy_rope(HY_RHOE,i1,i2,i3)
           
           hy_flux(HY_MASS:HY_ENER,i1,i2,i3) = &
           hy_flux(HY_MASS:HY_ENER,i1,i2,i3) &
           +cvisc*(/hy_rope(HY_DENS,i1-1,i2,i3)                 - hy_rope(HY_DENS,i1,i2,i3)&
           ,        hy_rope(HY_DENS,i1-1,i2,i3)*hy_rope(HY_VELX,i1-1,i2,i3) - hy_rope(HY_DENS,i1,i2,i3)*hy_rope(HY_VELX,i1,i2,i3)&
           ,        hy_rope(HY_DENS,i1-1,i2,i3)*hy_rope(HY_VELY,i1-1,i2,i3) - hy_rope(HY_DENS,i1,i2,i3)*hy_rope(HY_VELY,i1,i2,i3)&
           ,        hy_rope(HY_DENS,i1-1,i2,i3)*hy_rope(HY_VELZ,i1-1,i2,i3) - hy_rope(HY_DENS,i1,i2,i3)*hy_rope(HY_VELZ,i1,i2,i3)&
           ,        VenerLo                         - VenerHi/)
#ifdef SPARK_GLM
           hy_flux(HY_FMGX:HY_FPSI,i1,i2,i3) = &
           hy_flux(HY_FMGX:HY_FPSI,i1,i2,i3) &
           +cvisc*(/hy_rope(HY_MAGX,i1-1,i2,i3)                 - hy_rope(HY_MAGX,i1,i2,i3)&
           ,        hy_rope(HY_MAGY,i1-1,i2,i3)                 - hy_rope(HY_MAGY,i1,i2,i3)&
           ,        hy_rope(HY_MAGZ,i1-1,i2,i3)                 - hy_rope(HY_MAGZ,i1,i2,i3)&
           ,        hy_rope(HY_PSIB,i1-1,i2,i3)                 - hy_rope(HY_PSIB,i1,i2,i3)/)
#endif
           ! Here, we compute the species and mass scalar
           ! fluxes based on the density flux and the hy_reconstructed
           ! mass scalar interface values
           
           leftState => hy_uPlus
           rightState => hy_uMinus

#if NSPECIES+NMASS_SCALARS>0
           if (hy_flux(HY_MASS ,i1,i2,i3) > 0.) then
           hy_flux(HY_NUM_FLUX+1:NFLUXES ,i1,i2,i3) = leftState(HY_NUM_VARS+1:NRECON @M hy_m123)*hy_flux(HY_MASS ,i1,i2,i3)
           else
           hy_flux(HY_NUM_FLUX+1:NFLUXES ,i1,i2,i3) = rightState(HY_NUM_VARS+1:NRECON ,i1,i2,i3)*hy_flux(HY_MASS ,i1,i2,i3)
           end if
#endif
           nullify(leftState)
           nullify(rightState)
           end do
           end do
           end do
           ! ***************
           ! Fluxes computed for one face of this zone
           ! Save the fluxes
           ! ***************
           !$omp target teams distribute parallel do collapse(3) &
           !$omp private(i1,i2,i3,i_s,j_s,k_s) shared(dir,hy_dlim,gCells,hy_flx,hy_fly,hy_flz,hy_flux,klim) default(none)
           
           do i3=klim(LOW,KAXIS),klim(HIGH,KAXIS)
           do i2=klim(LOW,JAXIS),klim(HIGH,JAXIS)
           do i1=klim(LOW,IAXIS),klim(HIGH,IAXIS)      
           i_s = -1+i1+hy_dlim(LOW,1)-gCells(1)
           j_s = -1+i2+hy_dlim(LOW,2)-gCells(2)
           k_s = -1+i3+hy_dlim(LOW,3)-gCells(3)
           select case(dir)
           case(IAXIS)
           hy_flx(:,i_s,j_s,k_s) = hy_flux(:,i1,i2,i3)
           case (JAXIS)
           hy_fly(:,j_s,i_s,k_s) = hy_flux(:,i1,i2,i3)
           case (KAXIS)
           hy_flz(:,j_s,k_s,i_s) = hy_flux(:,i1,i2,i3)
           end select
           end do
           end do
           end do

           
        !release pointers
        
     
     
    
  end do ! dir
  
  !$omp end target data
 


contains

  !!Set loop dimensions based on the direction of the pencil set

  !! the 'pencil' holds a 1D array of the solution data to be operated
  !! on.  It is unrolled this way so that all data that are needed for
  !! interpolation and flux calculation are truly contiguous in memory
  !! space.  And so that it all fits in cache at once...
  !! The maximum amount of calculation is done on these data prior to
  !! reseting the pencil data to a new ray through the block.
  !! If I were talented at ASCII art, I would make a diagram...


  !~ Flattening has not been tested yet in FLASH5, only 1D & 2D runs so far.
  subroutine flattening(limits)
    !! This follows Miller & Colella 2002
    use Hydro_data, ONLY : hy_starState, hy_flat3d
    implicit none
    integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
    !real, intent(OUT) :: hy_flat3d(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    !real :: hy_flatTilde(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
    real :: hy_flatTilde(NDIM,limits(LOW,IAXIS):limits(HIGH,IAXIS),&
                              limits(LOW,JAXIS):limits(HIGH,JAXIS),&
                              limits(LOW,KAXIS):limits(HIGH,KAXIS))
    real :: beta, Z
    real, parameter :: betaMin = 0.75, betaMax = 0.85
    real, parameter :: Zmin = 0.25, Zmax = 0.75
    integer :: i,j,k, kx, ky, kz
    kx = 1
    ky = 0
    kz = 0
#if NDIM>1
    ky = 1
#if NDIM==3
    kz = 1
#endif
#endif
    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
             !1e-99 used to be TINY(1.0) but that produced Seg Faults
             beta = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
                  / max(1e-99, abs(hy_starState(PRES_VAR,i+2,j,k)-hy_starState(PRES_VAR,i-2,j,k)))
             Z    = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
             hy_flatTilde(IAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
             if (hy_starState(VELX_VAR,i+1,j,k)<hy_starState(VELX_VAR,i,j,k)) then
                hy_flatTilde(IAXIS,i,j,k) = max(hy_flatTilde(IAXIS,i,j,k), &
                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
             else
                hy_flatTilde(IAXIS,i,j,k) = 1.0
             end if
#if NDIM>1
             beta = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
                  / max(1.e-99,abs(hy_starState(PRES_VAR,i,j+2,k)-hy_starState(PRES_VAR,i,j-2,k)))
             Z    = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
             hy_flatTilde(JAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
             if (hy_starState(VELY_VAR,i,j+1,k)<hy_starState(VELY_VAR,i,j,k)) then
                hy_flatTilde(JAXIS,i,j,k) = max(hy_flatTilde(JAXIS,i,j,k), &
                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
             else
                hy_flatTilde(JAXIS,i,j,k) = 1.0
             end if
#if NDIM==3
             beta = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
                  / max(1.e-99,abs(hy_starState(PRES_VAR,i,j,k+2)-hy_starState(PRES_VAR,i,j,k-2)))
             Z    = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
             hy_flatTilde(KAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
             if (hy_starState(VELZ_VAR,i,j,k+1)<hy_starState(VELZ_VAR,i,j,k)) then
                hy_flatTilde(KAXIS,i,j,k) = max(hy_flatTilde(KAXIS,i,j,k), &
                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
             else
                hy_flatTilde(KAXIS,i,j,k) = 1.0
             end if
#endif
#endif
          end do
       end do
    end do
    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
             hy_flat3d(i,j,k) = minval(hy_flatTilde(1:NDIM,i,j,k))
#ifdef FLAT_VAR
             hy_starState(FLAT_VAR,i,j,k) = hy_flat3d(i,j,k)
#endif
          end do
       end do
    end do
  end subroutine flattening

end subroutine hy_rk_getFaceFlux

real function minmod(a,b)
  implicit none
  real :: a,b
  minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
end function minmod

