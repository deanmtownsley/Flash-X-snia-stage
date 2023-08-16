!!****if* source/physics/Hydro/HydroMain/unsplit/Hydro_computeDt
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
!!  Hydro_computeDt
!!
!!
!! DESCRIPTION
!!
!!  This routine computes the timestep limiter for the Unsplit Hydro solver.
!!  For more details see the documentation of the NULL implementation
!!
!!***

!!REORDER(4): U

Subroutine Hydro_computeDt( tileDesc,       &
                            x, dx, uxgrid, &
                            y, dy, uygrid, &
                            z, dz, uzgrid, &
                            blkLimits, blkLimitsGC, &
                            U,  dtCheck, dtMinLoc,  &
                            extraInfo)


#include "Simulation.h"
#include "constants.h"

  use Hydro_data, ONLY: hy_meshMe, hy_useHydro, hy_updateHydroFluxes, hy_hydroComputeDtOption, &
       hy_dtminValid, hy_hydroComputeDtFirstCall, hy_restart, &
       hy_dtmin, hy_dtminLoc, hy_dtminCfl, &
       hy_cfl, hy_cfl_original, hy_cflStencil, &
       hy_dref, hy_eref, hy_pref, hy_vref, hy_bref, &
       hy_geometry, hy_units, hy_useVaryingCFL
  use Driver_interface, ONLY : Driver_abort
  use Grid_tile, ONLY : Grid_tile_t 
  implicit none

  !! Arguments type declaration ------------------------------------------
  type(Grid_tile_t), intent(IN) :: tileDesc
  integer,dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC

  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)), intent(IN) :: x, dx, uxgrid
  real, dimension(blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)), intent(IN) :: y, dy, uygrid
  real, dimension(blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), intent(IN) :: z, dz, uzgrid


  real, pointer         :: U(:,:,:,:)
  real,   intent(INOUT) :: dtCheck
  integer,intent(INOUT) :: dtMinLoc(5)
  real, OPTIONAL,intent(INOUT) :: extraInfo
  !! ----------------------------------------------------------------------
  integer :: i, j, k, temploc(5)
  integer :: imS,ipS,jmS,jpS,kmS,kpS
  real    :: sndspd2, delxinv, delyinv, delzinv, dt_temp, dt_ltemp
  real    :: cfx2,cfy2,cfz2,bbx2,bby2,bbz2,b2
  real    :: localCfl, tempCfl, dtCflLoc
  integer, dimension(2,MDIM) :: bcs

  !! Case 1: we exit this routine if not needed.
  if ((.not. hy_useHydro) .or. (.not. hy_updateHydroFluxes)) return

  !! Case 2a: we pass the already computed dt information to Driver_computeDt
  !!          if it appears valid and the saved location's block ID matches the block argument.
  !!         In this case, hydro dt gets computed in either 
  !!          (i) hy_energyFix (hy_hydroComputeDtOption=0), or 
  !!         (ii) hy_getFaceFlux (hy_hydroComputeDtOption=1)

  if ((hy_hydroComputeDtOption .ne. -1) .and. &
       hy_dtminValid .and. &
      (.not. hy_hydroComputeDtFirstCall .OR. hy_restart)) then
     dtCflLoc = hy_cfl
     if ( hy_dtmin < dtCheck .AND. hy_dtminloc(4) == tileDesc%level) then
        dtCheck  = hy_dtmin
        dtMinLoc = hy_dtminloc(1:5)
        dtCflLoc = hy_dtminCfl
     endif

  !! Case 2: we simply pass the already computed dt information to Driver_computeDt.
  !!         In this case, hydro dt gets computed in either 
  !!          (i) hy_energyFix (hy_hydroComputeDtOption=0), or 
  !!         (ii) hy_getFaceFlux (hy_hydroComputeDtOption=1)
  else if ((hy_hydroComputeDtOption .ne. -1) .and. &
       hy_dtminValid .and. &
      (.not. hy_hydroComputeDtFirstCall)) then
     if ( hy_dtmin < dtCheck ) then
        dtCheck  = hy_dtmin
        dtMinLoc = hy_dtminloc(1:5)
     endif

     if (present(extraInfo)) then
        extraInfo = 0.
        if (hy_useVaryingCFL) then
           extraInfo = hy_cfl_original
           if (hy_cfl <= extraInfo) then
              extraInfo = hy_cfl
           endif
           dtCflLoc = extraInfo
        else
           extrainfo = 0.
        endif
     endif

  else
     !! Case 3: hy_hydroComputeDtOption = -1
     !!         We perform the global loop to compute hydro dt in the old way.
     !!         This implementation provides the full geometry supports for computing dt.
     !! NOTE: we always call this global compute dt for the very first step.
     if (hy_hydroComputeDtFirstCall) hy_hydroComputeDtFirstCall = .false.
     dt_temp    = 0.
     temploc(:) = 0

     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        !! Conversion for unitSystem if needed ---------------------------------------
        U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)/hy_dref
        U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)/hy_eref
        U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)/hy_pref
        U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)/hy_vref
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
        U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)/hy_bref
#endif
        !! ---------------------------------------------------------------------------
     endif


     dtCflLoc = hy_cfl_original
     tempCfl  = hy_cfl_original
#ifndef CFL_VAR
     localCfl = hy_cfl
#endif

     delxinv = 1.0/dx(blkLimits(LOW,IAXIS))
     if (NDIM > 1) &
     delyinv = 1.0/dy(blkLimits(LOW,JAXIS))
     if (NDIM > 2) &
     delzinv = 1.0/dz(blkLimits(LOW,KAXIS))
     if (hy_geometry == POLAR) & !Polar in 3D (that's a no no)
          call Driver_abort("[Hydro_computeDt] ERROR: Polar geometry not supported in 3D")

     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)

#ifdef BDRY_VAR /* Do not compute time step dt when in the solid boundary cells */
              if (U(BDRY_VAR,i,j,k) .LE. 0.0) then
#endif
#ifdef CFL_VAR
                 if (hy_cflStencil<1) then
                    localCfl = U(CFL_VAR,i,j,k)
                 else
                    call tileDesc%faceBCs(bcs)
                    imS=max(blkLimitsGC(LOW,IAXIS), i-hy_cflStencil)
                    ipS=min(blkLimitsGC(HIGH,IAXIS),i+hy_cflStencil)
                    
                    if (bcs(LOW,IAXIS) <= PARAMESH_PHYSICAL_BOUNDARY) &
                         imS = max(blkLimits(LOW,IAXIS), imS)
                    if (bcs(HIGH,IAXIS) <= PARAMESH_PHYSICAL_BOUNDARY) &
                         ipS = min(blkLimits(HIGH,IAXIS), ipS)
#if NDIM > 1
                    jmS=max(blkLimitsGC(LOW,JAXIS), j-hy_cflStencil)
                    jpS=min(blkLimitsGC(HIGH,JAXIS),j+hy_cflStencil)
                    
                    if (bcs(LOW,JAXIS) <= PARAMESH_PHYSICAL_BOUNDARY) &
                         jmS = max(blkLimits(LOW,JAXIS), jmS)
                    if (bcs(HIGH,JAXIS) <= PARAMESH_PHYSICAL_BOUNDARY) &
                         jpS = min(blkLimits(HIGH,JAXIS), jpS)
#else
                    jmS = 1; jpS=1
#endif
#if NDIM > 2
                    kmS=max(blkLimitsGC(LOW,KAXIS), k-hy_cflStencil)
                    kpS=min(blkLimitsGC(HIGH,KAXIS),k+hy_cflStencil)
                    
                    if (bcs(LOW,KAXIS) <= PARAMESH_PHYSICAL_BOUNDARY) &
                         kmS = max(blkLimits(LOW,KAXIS), kmS)
                    if (bcs(HIGH,KAXIS) <= PARAMESH_PHYSICAL_BOUNDARY) &
                         kpS = min(blkLimits(HIGH,KAXIS), kpS)
#else
                    kmS = 1; kpS=1
#endif
                    localCfl = minval(U(CFL_VAR,imS:ipS,jmS:jpS,kmS:kpS))
                 end if
#endif
                 sndspd2 = U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)
                 cfx2    = sndspd2
                 cfy2    = sndspd2
                 cfz2    = sndspd2

#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD) /*compute additional magneto-acoustic speeds for MHD */
                 bbx2 = U(MAGX_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 bby2 = U(MAGY_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 bbz2 = U(MAGZ_VAR,i,j,k)**2/U(DENS_VAR,i,j,k)
                 b2   = bbx2 + bby2 + bbz2
                 sndspd2= U(GAMC_VAR,i,j,k)*U(PRES_VAR,i,j,k)/U(DENS_VAR,i,j,k)

                 cfx2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bbx2))
                 cfy2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bby2))
                 cfz2 = .5*((sndspd2+b2)+sqrt((sndspd2+b2)**2-4.*sndspd2*bbz2))
#endif

                 ! For other geometry supports
                 if (hy_geometry == CYLINDRICAL) then
#if NDIM > 2
                    delzinv = 1.0/(x(i)*dz(k))           ! z is phi
#endif
                 elseif (hy_geometry == SPHERICAL) then
#if NDIM > 1
                    delyinv = 1.0/(x(i)*dy(j))           ! y is theta
#if NDIM > 2
                    delzinv = 1.0/(x(i)*sin(y(j))*dz(k)) ! z is phi
#endif
#endif
                 endif


                 dt_ltemp = (abs(U(VELX_VAR,i,j,k)-uxgrid(i))+sqrt(cfx2))*delxinv
                 if (NDIM > 1) dt_ltemp = max(dt_ltemp,(abs(U(VELY_VAR,i,j,k)-uygrid(j))+sqrt(cfy2))*delyinv)
                 if (NDIM > 2) dt_ltemp = max(dt_ltemp,(abs(U(VELZ_VAR,i,j,k)-uzgrid(k))+sqrt(cfz2))*delzinv)
                 if (dt_ltemp * tempCfl > dt_temp * localCfl) then
                    dt_temp    = dt_ltemp
                    tempCfl    = localCfl
                    temploc(1) = i
                    temploc(2) = j
                    temploc(3) = k
                    temploc(4) = tileDesc%level
                    temploc(5) = hy_meshMe
                 endif
#ifdef BDRY_VAR
              endif
#endif
           enddo
        enddo
     enddo


     if (dt_temp .NE. 0.0) then
        dt_temp = tempCfl / dt_temp
     else
        dt_temp = huge(1.0)
     end if
     if (dt_temp < dtCheck) then
        dtCheck = dt_temp
        dtMinLoc = temploc
        dtCflLoc = tempCfl
     endif

 
     if ( hy_units .NE. "NONE" .and. hy_units .NE. "none" ) then
        !! Conversion for unitSystem if needed ---------------------------------------
        U(DENS_VAR,:,:,:) = U(DENS_VAR,:,:,:)*hy_dref
        U(ENER_VAR,:,:,:) = U(ENER_VAR,:,:,:)*hy_eref
        U(PRES_VAR,:,:,:) = U(PRES_VAR,:,:,:)*hy_pref
        U(VELX_VAR:VELZ_VAR,:,:,:) = U(VELX_VAR:VELZ_VAR,:,:,:)*hy_vref
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
        U(MAGX_VAR:MAGZ_VAR,:,:,:) = U(MAGX_VAR:MAGZ_VAR,:,:,:)*hy_bref
#endif
        !! ---------------------------------------------------------------------------
     endif

  endif ! end if of if (hy_hydroComputeDtOption .ne. -1) then



  !! For the purpose of having screen output for varying CFL
  if (present(extraInfo)) then
     extraInfo = 0.
     if (hy_useVaryingCFL) then
        extraInfo = hy_cfl_original
        if (dtCflLoc <= extraInfo) then
           extraInfo = dtCflLoc
        endif
     else
        extrainfo = 0.
     endif
  endif

  if(dtCheck <= 0.0) then
     call Driver_abort("[Hydro]: Computed dt is not positive! Aborting!")
  endif
  return

End Subroutine Hydro_computeDt


