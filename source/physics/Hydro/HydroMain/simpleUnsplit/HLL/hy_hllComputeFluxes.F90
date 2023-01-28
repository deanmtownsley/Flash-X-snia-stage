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
! Note: the following arrays need to be spelled exactly like this in the code below,
!       preserving case.
!!REORDER(4): Uin, fl[XYZ], auxC

#ifdef DEBUG_ALL
#define DEBUG_UHD
#endif

#include "UHD.h"

Subroutine hy_hllComputeFluxes ( tileLimits, Uin, plo, flX, flY, flZ, loFl, del, dt )
  use Hydro_data,        ONLY : hy_fluxCorrect,      &
                                hy_useGravity,       &
                                hy_updateHydroFluxes
  use Driver_interface,  ONLY : Driver_abort
  use Logfile_interface, ONLY : Logfile_stampVarMask

  implicit none

  !! ---- Argument List ----------------------------------
  integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
  integer, intent(IN)  :: plo(*)
  real,    intent(IN)  :: UIN(plo(1):,plo(2):,plo(3):,plo(4):)  !CAPITALIZATION INTENTIONAL!
  integer, intent(IN)  :: loFl(*)
  real,    intent(OUT) :: FLX(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: FLY(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(OUT) :: FLZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):) !CAPITALIZATION INTENTIONAL!
  real,    intent(IN)  :: del(1:MDIM)
  real,    intent(IN)  :: dt
  !! -----------------------------------------------------

#ifdef INDEXREORDER
  integer, parameter :: iX = 1
#else
  integer, parameter :: iX = 2
#endif
  integer :: i,j,k
  integer :: is,js,ks
  integer :: iLastX, iLastY, iLastZ
  integer :: iL,iR, jL, jR, kL, kR
  real :: dtdx, dtdy, dtdz, vn
  real :: c, sL, sR
  real :: sRsL, vL, vR

  real, pointer, dimension(:,:,:,:) :: auxC

  nullify(auxC)

  !! End of data declaration ***********************************************
#ifdef DEBUG_UHD
98 format(A4,'(',I3,':   ,',   I3,':   ,',   I3,':   ,',   I3,':   )')
99 format(A4,'(',I3,':',I3,',',I3,':',I3,',',I3,':',I3,',',I3,':',I3,')')
  print *, "plo" ,plo(1:MDIM+1)
  print 98,"Uin" ,(plo(i),i=1,4)
  print 99,"Uin" ,(lbound(Uin ,i),ubound(Uin ,i),i=1,4)
  print 99,"Uout",(lbound(Uout,i),ubound(Uout,i),i=1,4)
  print*,'tileLim:',tileLimits
#endif

  if (hy_fluxCorrect) then
     call Driver_abort("hy_hllUnsplit: flux correction is not implemented!")
  end if

  if (hy_useGravity) then
     call Driver_abort("hy_hllUnsplit: support for gravity not implemented!")
  end if

  if (.NOT.hy_updateHydroFluxes) then
     return
  end if

  dtdx = dt / del(IAXIS)
  if (NDIM > 1) dtdy = dt / del(JAXIS)
  if (NDIM > 2) dtdz = dt / del(KAXIS)

  ! For each relevant direction We set indicators to determine whether
  ! fluxes on the last face (to the right of the last cell in the
  ! cell-centered range given in tileLimits) need to be computed or
  ! not. They need to be computed only if the face is on the boundary
  ! (between the last interior in the first guard cell). If proper
  ! tiling is in effect, we want to avoid computing the same flux
  ! twice if the last face of some tile coincides with the first face
  ! of the next tile.
  ! The following logic assumes that Uin is sized such that the high
  ! bound of its relevant indices is exactly NGUARD cells to the
  ! right of the last interior cell.
  iLastX = 0; iLastY = 0; iLastZ = 0
  if (tileLimits(HIGH,IAXIS) == ubound(Uin,iX  )-NGUARD) iLastX = 1
#if NDIM > 1
  if (tileLimits(HIGH,JAXIS) == ubound(Uin,iX+1)-NGUARD) iLastY = K2D
#endif
#if NDIM > 2
  if (tileLimits(HIGH,KAXIS) == ubound(Uin,iX+2)-NGUARD) iLastZ = K3D
#endif

  allocate(auxC(1,tileLimits(LOW,IAXIS)-1  :tileLimits(HIGH,IAXIS)+iLastX, &
                  tileLimits(LOW,JAXIS)-K2D:tileLimits(HIGH,JAXIS)+iLastY, &
                  tileLimits(LOW,KAXIS)-K3D:tileLimits(HIGH,KAXIS)+iLastX) )

  !! ************************************************************************
  !! Calculate Riemann (interface) states
  !  No equivalent really to  call hy_uhd_getRiemannState(blockID,blkLimits,blkLimitsGC,dt,del)

  !! calculate sound speed
  do k = tileLimits(LOW,KAXIS)-K3D,tileLimits(HIGH,KAXIS)+iLastZ
     do j = tileLimits(LOW,JAXIS)-K2D,tileLimits(HIGH,JAXIS)+iLastY
        do i = tileLimits(LOW,IAXIS)-1,tileLimits(HIGH,IAXIS)+iLastX
           c = sqrt(Uin(GAMC_VAR,i,j,k)*Uin(PRES_VAR,i,j,k)/Uin(DENS_VAR,i,j,k))
           auxC(1, i,j,k) = c
        end do
     end do
  end do

  !! ************************************************************************
  !! Calculate Godunov fluxes

  !  instead of  call hy_uhd_getFaceFlux(blockID,blkLimits,blkLimitsGC,datasize,del, ...)

  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)+iLastX
           sL = min(Uin(VELX_VAR,i-1,j,k)-auxC(1, i-1,j,k), Uin(VELX_VAR,i,j,k)-auxC(1, i,j,k))
           sR = max(Uin(VELX_VAR,i-1,j,k)+auxC(1, i-1,j,k), Uin(VELX_VAR,i,j,k)+auxC(1, i,j,k))
           sRsL = sR - sL
           if (sL > 0.0) then
              vn = Uin(VELX_VAR,i-1,j,k)
              is = i-1
              iL = i-1; iR=i-1
           else if (sR < 0.0) then
              vn = Uin(VELX_VAR,i,j,k)
              is = i
              iL = i; iR=i
           else
              vn = (Uin(VELX_VAR,i-1,j,k)+Uin(VELX_VAR,i,j,k)) * 0.5
              is = i
              iL = i-1; iR=i
              if (vn>0.0) is = is-1
           end if
           vL = Uin(VELX_VAR,iL,j,k);  vR = Uin(VELX_VAR,iR,j,k)
           if (iL==iR) then
              flX(HY_DENS_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k)
              flX(HY_XMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(VELX_VAR,is,j,k)
              flX(HY_XMOM_FLUX,i,j,k) = flX(HY_XMOM_FLUX,i,j,k) &
                   + Uin(PRES_VAR,is,j,k)
              flX(HY_YMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(VELY_VAR,is,j,k)
              flX(HY_ZMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(VELZ_VAR,is,j,k)
              flX(HY_ENER_FLUX,i,j,k) = vn * Uin(DENS_VAR,is,j,k) * Uin(ENER_VAR,is,j,k)
              flX(HY_ENER_FLUX,i,j,k) = flX(HY_ENER_FLUX,i,j,k) &
                   + vn * Uin(PRES_VAR,is,j,k)
           else
              flX(HY_DENS_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k) - sL * vR * Uin(DENS_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k))) / sRsL
              flX(HY_XMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(VELX_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(VELX_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(VELX_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELX_VAR,iL,j,k)) )/sRsL
              flX(HY_XMOM_FLUX,i,j,k) = flX(HY_XMOM_FLUX,i,j,k) &
                   + (sR * Uin(PRES_VAR,iL,j,k) - sL * Uin(PRES_VAR,iR,j,k)) / sRsL
              flX(HY_YMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(VELY_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(VELY_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(VELY_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELY_VAR,iL,j,k)) )/sRsL
              flX(HY_ZMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(VELZ_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(VELZ_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(VELZ_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(VELZ_VAR,iL,j,k)) )/sRsL
              flX(HY_ENER_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,iL,j,k)*Uin(ENER_VAR,iL,j,k) &
                   - sL * vR * Uin(DENS_VAR,iR,j,k)*Uin(ENER_VAR,iR,j,k) &
                   + sR*sL*(Uin(DENS_VAR,iR,j,k)*Uin(ENER_VAR,iR,j,k) - Uin(DENS_VAR,iL,j,k)*Uin(ENER_VAR,iL,j,k)))/sRsL
              flX(HY_ENER_FLUX,i,j,k) = flX(HY_ENER_FLUX,i,j,k) &
                   + (sR * vL * Uin(PRES_VAR,iL,j,k) - sL * vR * Uin(PRES_VAR,iR,j,k)) / sRsL
           end if
           flX(HY_DENS_FLUX,i,j,k) = flX(HY_DENS_FLUX,i,j,k) * dtdx
           flX(HY_XMOM_FLUX,i,j,k) = flX(HY_XMOM_FLUX,i,j,k) * dtdx
           flX(HY_YMOM_FLUX,i,j,k) = flX(HY_YMOM_FLUX,i,j,k) * dtdx
           flX(HY_ZMOM_FLUX,i,j,k) = flX(HY_ZMOM_FLUX,i,j,k) * dtdx
           flX(HY_ENER_FLUX,i,j,k) = flX(HY_ENER_FLUX,i,j,k) * dtdx
        end do
     end do
  end do
#if NDIM > 1
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)+iLastY
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           sL = min(Uin(VELY_VAR,i,j-1,k)-auxC(1, i,j-1,k), Uin(VELY_VAR,i,j,k)-auxC(1, i,j,k))
           sR = max(Uin(VELY_VAR,i,j-1,k)+auxC(1, i,j-1,k), Uin(VELY_VAR,i,j,k)+auxC(1, i,j,k))
           sRsL = sR - sL
           if (sL > 0.0) then
              vn = Uin(VELY_VAR,i,j-1,k)
              js = j-1
              jL = j-1; jR=j-1
           else if (sR < 0.0) then
              vn = Uin(VELY_VAR,i,j,k)
              js = j
              jL = j; jR=j
           else
              vn = (Uin(VELY_VAR,i,j-1,k)+Uin(VELY_VAR,i,j,k)) * 0.5
              js = j
              jL = j-1; jR=j
              if (vn>0.0) js = js-1
           end if
           vL = Uin(VELY_VAR,i,jL,k);  vR = Uin(VELY_VAR,i,jR,k)
           if (jL==jR) then
              flY(HY_DENS_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k)
              flY(HY_XMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(VELX_VAR,i,js,k)
              flY(HY_YMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(VELY_VAR,i,js,k)
              flY(HY_YMOM_FLUX,i,j,k) = flY(HY_YMOM_FLUX,i,j,k) &
                   + Uin(PRES_VAR,i,js,k)
              flY(HY_ZMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(VELZ_VAR,i,js,k)
              flY(HY_ENER_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,js,k) * Uin(ENER_VAR,i,js,k)
              flY(HY_ENER_FLUX,i,j,k) = flY(HY_ENER_FLUX,i,j,k) &
                   + vn * Uin(PRES_VAR,i,js,k)
           else
              flY(HY_DENS_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k) - sL * vR * Uin(DENS_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k))) / sRsL
              flY(HY_XMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(VELX_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(VELX_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(VELX_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELX_VAR,i,jL,k)) )/sRsL
              flY(HY_YMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(VELY_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(VELY_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(VELY_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELY_VAR,i,jL,k)) )/sRsL
              flY(HY_YMOM_FLUX,i,j,k) = flY(HY_YMOM_FLUX,i,j,k) &
                   + (sR * Uin(PRES_VAR,i,jL,k) - sL * Uin(PRES_VAR,i,jR,k)) / sRsL
              flY(HY_ZMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(VELZ_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(VELZ_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(VELZ_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(VELZ_VAR,i,jL,k)) )/sRsL
              flY(HY_ENER_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,jL,k)*Uin(ENER_VAR,i,jL,k) &
                   - sL * vR * Uin(DENS_VAR,i,jR,k)*Uin(ENER_VAR,i,jR,k) &
                   + sR*sL*(Uin(DENS_VAR,i,jR,k)*Uin(ENER_VAR,i,jR,k) - Uin(DENS_VAR,i,jL,k)*Uin(ENER_VAR,i,jL,k)))/sRsL
              flY(HY_ENER_FLUX,i,j,k) = flY(HY_ENER_FLUX,i,j,k) &
                   + (sR * vL * Uin(PRES_VAR,i,jL,k) - sL * vR * Uin(PRES_VAR,i,jR,k)) / sRsL
           end if
           flY(HY_DENS_FLUX,i,j,k) = flY(HY_DENS_FLUX,i,j,k) * dtdy
           flY(HY_XMOM_FLUX,i,j,k) = flY(HY_XMOM_FLUX,i,j,k) * dtdy
           flY(HY_YMOM_FLUX,i,j,k) = flY(HY_YMOM_FLUX,i,j,k) * dtdy
           flY(HY_ZMOM_FLUX,i,j,k) = flY(HY_ZMOM_FLUX,i,j,k) * dtdy
           flY(HY_ENER_FLUX,i,j,k) = flY(HY_ENER_FLUX,i,j,k) * dtdy
        end do
     end do
  end do
#else
  flY(:,:,:,:) = 0.0      !avoid compiler warning or error for intent(OUT) dummy
#endif
#if NDIM > 2
  do k = tileLimits(LOW,KAXIS),tileLimits(HIGH,KAXIS)+iLastZ
     do j = tileLimits(LOW,JAXIS),tileLimits(HIGH,JAXIS)
        do i = tileLimits(LOW,IAXIS),tileLimits(HIGH,IAXIS)
           sL = min(Uin(VELZ_VAR,i,j,k-1)-auxC(1, i,j,k-1), Uin(VELZ_VAR,i,j,k)-auxC(1, i,j,k))
           sR = max(Uin(VELZ_VAR,i,j,k-1)+auxC(1, i,j,k-1), Uin(VELZ_VAR,i,j,k)+auxC(1, i,j,k))
           sRsL = sR - sL
           if (sL > 0.0) then
              vn = Uin(VELZ_VAR,i,j,k-1)
              ks = k-1
              kL = k-1; kR = k-1
           else if (sR < 0.0) then
              vn = Uin(VELZ_VAR,i,j,k)
              ks = k
              kL = k; kR = k
           else
              vn = (Uin(VELZ_VAR,i,j,k-1)+Uin(VELZ_VAR,i,j,k)) * 0.5
              ks = k
              kL = k-1; kR = k
              if (vn>0.0) ks = ks-1
           end if
           vL = Uin(VELZ_VAR,i,j,kL);  vR = Uin(VELZ_VAR,i,j,kR)
           if (kL==kR) then
              flZ(HY_DENS_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks)
              flZ(HY_XMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(VELX_VAR,i,j,ks)
              flZ(HY_YMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(VELY_VAR,i,j,ks)
              flZ(HY_ZMOM_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(VELZ_VAR,i,j,ks)
              flZ(HY_ZMOM_FLUX,i,j,k) = flZ(HY_ZMOM_FLUX,i,j,k) &
                   + Uin(PRES_VAR,i,j,ks)
              flZ(HY_ENER_FLUX,i,j,k) = vn * Uin(DENS_VAR,i,j,ks) * Uin(ENER_VAR,i,j,ks)
              flZ(HY_ENER_FLUX,i,j,k) = flZ(HY_ENER_FLUX,i,j,k) &
                   + vn * Uin(PRES_VAR,i,j,ks)
           else

              flZ(HY_DENS_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL) - sL * vR * Uin(DENS_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL))) / sRsL
              flZ(HY_XMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(VELX_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(VELX_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(VELX_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELX_VAR,i,j,kL)) )/sRsL
              flZ(HY_YMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(VELY_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(VELY_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(VELY_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELY_VAR,i,j,kL)) )/sRsL
              flZ(HY_ZMOM_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(VELZ_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(VELZ_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(VELZ_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(VELZ_VAR,i,j,kL)) )/sRsL
              flZ(HY_ZMOM_FLUX,i,j,k) = flZ(HY_ZMOM_FLUX,i,j,k) &
                   + (sR * Uin(PRES_VAR,i,j,kL) - sL * Uin(PRES_VAR,i,j,kR)) / sRsL
              flZ(HY_ENER_FLUX,i,j,k) = (  sR * vL * Uin(DENS_VAR,i,j,kL)*Uin(ENER_VAR,i,j,kL) &
                   - sL * vR * Uin(DENS_VAR,i,j,kR)*Uin(ENER_VAR,i,j,kR) &
                   + sR*sL*(Uin(DENS_VAR,i,j,kR)*Uin(ENER_VAR,i,j,kR) - Uin(DENS_VAR,i,j,kL)*Uin(ENER_VAR,i,j,kL)))/sRsL
              flZ(HY_ENER_FLUX,i,j,k) = flZ(HY_ENER_FLUX,i,j,k) &
                   + (sR * vL * Uin(PRES_VAR,i,j,kL) - sL * vR * Uin(PRES_VAR,i,j,kR)) / sRsL
           end if
           flZ(HY_DENS_FLUX,i,j,k) = flZ(HY_DENS_FLUX,i,j,k) * dtdz
           flZ(HY_XMOM_FLUX,i,j,k) = flZ(HY_XMOM_FLUX,i,j,k) * dtdz
           flZ(HY_YMOM_FLUX,i,j,k) = flZ(HY_YMOM_FLUX,i,j,k) * dtdz
           flZ(HY_ZMOM_FLUX,i,j,k) = flZ(HY_ZMOM_FLUX,i,j,k) * dtdz
           flZ(HY_ENER_FLUX,i,j,k) = flZ(HY_ENER_FLUX,i,j,k) * dtdz
        end do
     end do
  end do
#else
  flZ(:,:,:,:) = 0.0      !avoid compiler warning or error for intent(OUT) dummy
#endif

  deallocate(auxC)

End Subroutine hy_hllComputeFluxes

