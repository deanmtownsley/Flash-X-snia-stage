!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_correctFluxes
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
!!  hy_rk_correctFluxes
!!
!!  SYNOPSIS
!!
!!  call hy_rk_correctFluxes (real, pointer :: Uin(:,:,:,:),
!!                            integer (IN)  :: blkLimits(:,:),
!!                            integer (IN)  :: blkLimitsGC(:,:),
!!                            integer (IN)  :: level,
!!                            real (IN)     :: hy_del(:),
!!                            real (IN)     :: dt)
!!     
!!
!!  DESCRIPTION
!!  Apply flux correction at fine/coarse boundaries. 
!!  The proper 'flux hy_del'--consistent with Section 3 of Berger&Colella(1989)--
!!  should be loaded into hy_fluxBuf[XYZ].  Because this algorithm applies flux 
!!  correction at every block interface (not just fine/coarse interfaces) the fluxes 
!!  NOT on fine/coarse interfaces must be zerod out. 
!!
!!
!!  ARGUMENTS
!!    Uin -- pointer to solution data
!!    blkLimits, blkLimitsGC -- index limits for interior and exterior of the tile
!!    level  -- the refine level of the block
!!    hy_del  --- dx, dy, dz
!!    dt - time step
!!
!!
!!***
!!Reorder(4):p_fluxBuf[XYZ],Uin
subroutine hy_rk_correctFluxes(Uin,blkLimits,BlklimitsGC,level,hy_del, dt)

  use Hydro_data, ONLY : hy_threadWithinBlock, hy_starState, &
       hy_smallE, hy_smalldens, hy_geometry, &
       hy_grav, hy_4piGinv, hy_alphaGLM, hy_C_hyp, hy_fluxCorVars, &
       hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ,hy_farea,hy_cvol,&
       hy_xCenter,hy_xLeft,hy_xRight,hy_eosData, hy_mfrac
  use Driver_interface, ONLY : Driver_abort
  use Eos_interface, ONLY : Eos_wrapped,Eos_getData,Eos_putData,Eos


  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"

  real, pointer :: Uin(:,:,:,:)
  integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  integer,intent(IN) :: level
  real,dimension(MDIM) :: hy_del
  real, intent(IN) :: dt

  integer, dimension(MDIM) :: lo, hi
  
  integer :: i,j,k,n,g


  real, pointer :: p_fluxBufX(:,:,:,:),p_fluxBufY(:,:,:,:),p_fluxBufZ(:,:,:,:)
  real, pointer :: Vstar(:)
  real :: dx, dy, dz
  real :: dFlux(NFLUXES)

  real :: eint, ekin, emag
  ! Geometry factors
  real :: facM, facP
  real :: dhdt, fac
 
  ! For EOS call
  !integer,dimension(MDIM) :: pos
  integer :: vecLen = 1 !b/c updated cell by cell
  integer, dimension(LOW:HIGH,MDIM)  :: range

  lo(:) = blkLimits(LOW,:)
  hi(:) = blkLimits(HIGH,:)

 
 
   dhdt = minval(hy_del(1:NDIM))/dt

  !~ hy_fluxBuf[XYZ] represents (sum(F_fine) - F_coarse) on 
  !~ coarse side of f/c boundary, 0 elsewhere

  nullify(p_fluxBufX);nullify(p_fluxBufY);nullify(p_fluxBufZ)
  !These pointers allow us to use hy_fluxBuf[XYZ] (whose limits
  !are hard coded in FBS mode) within the varying loop limits below

  p_fluxBufX(1:NFLUXES,lo(IAXIS):hi(IAXIS)+1,lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) => hy_fluxBufX
#if NDIM>1  
  p_fluxBufY(1:NFLUXES,lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS)+1,lo(KAXIS):hi(KAXIS)) => hy_fluxBufY
#if NDIM==3
  p_fluxBufZ(1:NFLUXES,lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)+1) => hy_fluxBufZ
#endif
#endif

 
!  !$omp parallel if (.FALSE.) &
!  !$omp default(none) &
!  !$omp firstprivate(vecLen)&
!  !$omp private(n,i,j,k,Vstar,facM,facP,range,hy_eosData,hy_mfrac,&
!  !$omp         emag,ekin,dx,dy,dz,fac,dFlux)&
!  !$omp shared(dt,Uin,hy_starState,p_fluxBufX,p_fluxBufY,p_fluxBufZ,hy_grav,&
!  !$omp        hy_fluxBufX,hy_fluxBufY,hy_fluxBufZ,hy_xCenter,hy_xLeft,hy_xRight,&
!  !$omp        blkLimits,blkLimitsGC,hy_alphaGLM, hy_C_hyp,&
!  !$omp        dhdt, hy_smalldens, hy_smallE,hy_del)

  ! Correct IAXIS sides
  !$omp do schedule(static) collapse(2)
  !Change limits to grownTile?
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        !! LOW side
        i = blkLimits(LOW,IAXIS)
        ! Point to old/intermediate states
        Vstar => Uin(:,i,j,k)
        ! Point to the correct fluxes
        !~ (-) sign b/c fluxes are leaving on the low side
        dFlux = -p_fluxBufX(:,i  ,j  ,k  )
        dx = hy_del(IAXIS)
        ! Get geometric factors and sources
        call geoFacs(i,j,k,facM,facP)
        fac = facM
        ! if (dFlux(HY_ENER) /= 0.) print *, 'b', Vstar(TEMP_VAR), dt/dx*fac*dFlux(HY_ENER)/(Vstar(ENER_VAR)*Vstar(DENS_VAR)), Vstar(VELX_VAR)
        !Update primitives (Vstar)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i
        range(LOW:HIGH,JAXIS) = j
        range(LOW:HIGH,KAXIS) = k
        !call Eos_getData(range,vecLen,tempData,CENTER,hy_eosData,hy_mfrac)
        call Eos_getData(range,vecLen,Uin,CENTER,hy_eosData,hy_mfrac)
        !call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        !call Eos_putData(range,vecLen,tempData,CENTER,hy_eosData)
        call Eos_putData(range,vecLen,Uin,CENTER,hy_eosData)
!!$        call Eos_wrapped(MODE_DENS_EI,range,Uin)
        ! if (dFlux(HY_ENER) /= 0.) print *, 'a', Vstar(TEMP_VAR), dt/dx*fac*dFlux(HY_ENER)/(Vstar(ENER_VAR)*Vstar(DENS_VAR)), Vstar(VELX_VAR)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side
        i = blkLimits(HIGH,IAXIS)
        ! Point to old/intermediate states
        Vstar => Uin(:,i,j,k)
        ! Point to the correct fluxes
        ! Positive b/c fluxes are entering on the high side
        dFlux = p_fluxBufX(:,i+1  ,j  ,k  )
        dx = hy_del(IAXIS)
        ! Get geometric factors and sources
        call geoFacs(i,j,k,facM,facP)
        fac = facP
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i 
        range(LOW:HIGH,JAXIS) = j 
        range(LOW:HIGH,KAXIS) = k
       
        call Eos_getData(range,vecLen,Uin,CENTER,hy_eosData,hy_mfrac)
        call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        call Eos_putData(range,vecLen,Uin,CENTER,hy_eosData)

        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#if NDIM>1
  ! Correct JAXIS sides
  !$omp do schedule(static) collapse(2)
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        !! LOW side
        j = blkLimits(LOW,JAXIS)
        ! Point to old/intermediate states
        Vstar => Uin(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -p_fluxBufY(:,i  ,j  ,k  )
        fac = 1.0
        dx = hy_del(JAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i 
        range(LOW:HIGH,JAXIS) = j 
        range(LOW:HIGH,KAXIS) = k

        call Eos_getData(range,vecLen,Uin,CENTER,hy_eosData,hy_mfrac)
        call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        call Eos_putData(range,vecLen,Uin,CENTER,hy_eosData)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side

        j = blkLimits(HIGH,JAXIS)
        ! Point to old/intermediate states
        Vstar => Uin(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = p_fluxBufY(:,i  ,j+1  ,k  )
        fac = 1.0
        dx = hy_del(JAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i 
        range(LOW:HIGH,JAXIS) = j 
        range(LOW:HIGH,KAXIS) = k
       
        call Eos_getData(range,vecLen,Uin,CENTER,hy_eosData,hy_mfrac)
        call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        call Eos_putData(range,vecLen,Uin,CENTER,hy_eosData)
        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#if NDIM==3
  ! Correct KAXIS sides
  !$omp do schedule(static) collapse(2)
  do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
     do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
        !! LOW side
        k = blkLimits(LOW,KAXIS)
        ! Point to old/intermediate states
        Vstar => Uin(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -p_fluxBufZ(:,i  ,j  ,k  )
        fac = 1.0
        dx = hy_del(KAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i
        range(LOW:HIGH,JAXIS) = j
        range(LOW:HIGH,KAXIS) = k
        
        call Eos_getData(range,vecLen,Uin,CENTER,hy_eosData,hy_mfrac)
        call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        call Eos_putData(range,vecLen,Uin,CENTER,hy_eosData)
        
        !release pointers
        nullify(Vstar)

        !! HIGH side
        k = blkLimits(HIGH,KAXIS)
        ! Point to old/intermediate states
        Vstar => Uin(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = p_fluxBufZ(:,i  ,j  ,k+1  )
        fac = 1.0
        dx = hy_del(KAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i
        range(LOW:HIGH,JAXIS) = j
        range(LOW:HIGH,KAXIS) = k
       
        call Eos_getData(range,vecLen,Uin,CENTER,hy_eosData,hy_mfrac)
        call Eos(MODE_DENS_EI,vecLen,hy_eosData,hy_mfrac)
        call Eos_putData(range,vecLen,Uin,CENTER,hy_eosData) 
        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#endif
#endif

!  !$omp end parallel

  nullify(p_fluxBufX);nullify(p_fluxBufY);nullify(p_fluxBufZ)


contains

  !!Update cell state based on flux hy_del
  subroutine correctZone(Vstar,dFlux,dt,dx,fac)
    implicit none
    real, pointer, intent(INOUT) :: Vstar(:)
    real, intent(IN) :: dFlux(NFLUXES), dt, dx, fac
    real :: Ustar(NFLUXES)

    ! Construct vectors of conserved variables
    Ustar(HY_MASS)         = Vstar(DENS_VAR)
    Ustar(HY_XMOM:HY_ZMOM) = Vstar(DENS_VAR)*Vstar(VELX_VAR:VELZ_VAR)
    Ustar(HY_ENER)         = Vstar(DENS_VAR)*Vstar(ENER_VAR)
    Ustar(HY_NUM_FLUX+1:NFLUXES) = Vstar(SPECIES_BEGIN:MASS_SCALARS_END)*Vstar(DENS_VAR)
#ifdef SPARK_GLM
    Ustar(HY_FMGX:HY_FMGZ) = Vstar(MAGX_VAR:MAGZ_VAR)
    Ustar(HY_ENER) = Ustar(HY_ENER)+0.5*dot_product(Vstar(MAGX_VAR:MAGZ_VAR),&
         Vstar(MAGX_VAR:MAGZ_VAR))  ! * density ??? [KC]
    Ustar(HY_FPSI) = Vstar(PSIB_VAR)
#endif

    ! Now correct conserved vector with flux hy_del
    ! The facP, facM definition is cracked. Need to know what side we are on
    Ustar = Ustar -dt/dx*(fac*dFlux)

    ! Update primitive variables
    emag = 0.0
#ifdef SPARK_GLM
    Vstar(MAGX_VAR:MAGZ_VAR) = Ustar(HY_FMGX:HY_FMGZ)
    ! Parabolic damping of PSI is applied to flux correction difference above
    Vstar(PSIB_VAR) = Ustar(HY_FPSI)
    emag = 0.5*dot_product(Vstar(MAGX_VAR:MAGZ_VAR),Vstar(MAGX_VAR:MAGZ_VAR))
    Vstar(MAGP_VAR) = emag
    Ustar(HY_ENER) = Ustar(HY_ENER) - emag
#endif
    Vstar(DENS_VAR)          = max(Ustar(HY_MASS),hy_smalldens)
    Vstar(VELX_VAR:VELZ_VAR) = Ustar(HY_XMOM:HY_ZMOM)/Vstar(DENS_VAR)
    Vstar(ENER_VAR)          = max(hy_smallE,Ustar(HY_ENER)/Vstar(DENS_VAR))

    ekin = .5*dot_product(Vstar(VELX_VAR:VELZ_VAR),Vstar(VELX_VAR:VELZ_VAR))
    Vstar(EINT_VAR) = max(hy_smallE,Vstar(ENER_VAR)-ekin)

    ! Divide partial densities by new mass densities to finalize
    ! update of new mass fractions.
    Vstar(SPECIES_BEGIN:MASS_SCALARS_END) = Ustar(HY_NUM_FLUX+1:NFLUXES)/Vstar(DENS_VAR)
  end subroutine correctZone

  !!Geometric factors for non-Cartesian geometries
  subroutine  geoFacs(i,j,k,facM,facP)
    implicit none
    integer, intent(IN) :: i,j,k
    real, intent(OUT) :: facM, facP
    real    :: alpha, dx_sph

    if (hy_geometry == CARTESIAN) then
       facM = 1.0; facP = 1.0
       return
    endif
    facM = hy_farea(i  ,j,k)*dx/hy_cvol(i,j,k)
    facP = hy_farea(i+1,j,k)*dx/hy_cvol(i,j,k)

    if (hy_xCenter(i) < 0.0) then
       facM = 0.
       facP = 0.
    end if
  end subroutine geoFacs

end subroutine hy_rk_correctFluxes
