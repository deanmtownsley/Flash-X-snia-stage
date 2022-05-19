!!****if* source/physics/Hydro/HydroMain/Spark/hy_rk_correctFluxes
!!
!!  NAME
!!
!!  hy_rk_correctFluxes
!!
!!  SYNOPSIS
!!
!!  call hy_rk_correctFluxes ( type(Grid_tile_t) :: blockDesc,
!!                              real, intent(IN) :: dt )
!!
!!  DESCRIPTION
!!  Apply flux correction at fine/coarse boundaries. 
!!  The proper 'flux deltas'--consistent with Section 3 of Berger&Colella(1989)--
!!  should be loaded into hy_fluxBuf[XYZ].  Because this algorithm applies flux 
!!  correction at every block interface (not just fine/coarse interfaces) the fluxes 
!!  NOT on fine/coarse interfaces must be zerod out. 
!!
!!
!!  ARGUMENTS
!!  blockDesc - block/tile descriptor
!!  dt - time step
!!***
!!Reorder(4):p_fluxBuf[XYZ],solnData
subroutine hy_rk_correctFluxes(blockDesc, dt)

  use Hydro_data, ONLY : hy_threadWithinBlock, hy_starState, &
       hy_smallE, hy_smalldens, hy_geometry, &
       hy_grav, hy_4piGinv, hy_alphaGLM, hy_C_hyp, hy_fluxCorVars, &
       hya_flxbx, hya_flxby, hya_flxbz
  use Driver_interface, ONLY : Driver_abort
  use Grid_interface, ONLY : Grid_getCellCoords, Grid_getCellFaceAreas, &
                             Grid_getCellVolumes
  use Eos_interface, ONLY : Eos_putData, Eos_getData, Eos
  use Grid_tile, ONLY : Grid_tile_t


  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"

  type(Grid_tile_t) :: blockDesc
  real, intent(IN) :: dt

  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
  integer, dimension(MDIM) :: lo, hi, loGC, hiGC
  integer :: xLoGC,yLoGC,zLoGC,xHiGC,yHiGC,zHiGC

  real, allocatable, dimension(:) :: xCenter, xLeft, xRight
  integer :: i,j,k,n,g, lev

  real, pointer :: solnData(:,:,:,:)
  real, pointer :: p_fluxBufX(:,:,:,:),p_fluxBufY(:,:,:,:),p_fluxBufZ(:,:,:,:)
  real, pointer :: Vstar(:)
  real :: dx, dy, dz, del(MDIM)
  real :: dFlux(NFLUXES)

  real :: eint, ekin, emag
  ! Geometry factors
  real :: facM, facP
  integer :: isize, jsize, ksize
  integer, dimension(MDIM) :: datasize
  real :: dhdt, fac
  real, allocatable, dimension(:,:,:) :: faceAreas, cellVolumes
 
  ! For EOS call
  !integer,dimension(MDIM) :: pos
  integer :: vecLen = 1 !b/c updated cell by cell
  integer, dimension(LOW:HIGH,MDIM)  :: range
  !Note this allocation works around the inability of FLASH5 to produce MAXCELLS constant in 
  !Simulation.h using preprocessors alone.  This may want to be changed for future.
  !For example the blkLimitsGC could be passed in as arguments.  This would reduce the 
  !amount of dynamic allocation but clutter the code.
  !real, dimension(NSPECIES*MAXCELLS) :: massFraction
  !real, dimension(EOS_NUM*MAXCELLS) :: eosData
  real, allocatable :: massFraction(:), eosData(:)

  blkLimits(:,:) = blockDesc%limits

  lo(:) = blkLimits(LOW,:)
  hi(:) = blkLimits(HIGH,:)

  blkLimitsGC(:,:) = blockDesc%blkLimitsGC
  loGC(:) = blkLimitsGC(LOW,:)
  hiGC(:) = blkLimitsGC(HIGH,:)
  !convenience indices
  xLoGC = loGC(IAXIS); xHiGC = hiGC(IAXIS)
  yLoGC = loGC(JAXIS); yHiGC = hiGC(JAXIS)
  zLoGC = loGC(KAXIS); zHiGC = hiGC(KAXIS)

  iSize = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
  jSize = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
  kSize = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
  datasize(1:MDIM) = blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
 
  allocate(massFraction(NSPECIES*MAXVAL(datasize)))
  allocate(eosData(EOS_NUM*MAXVAL(datasize)))
 
  nullify(solnData)
  
  call blockDesc%getDataPtr(solnData,CENTER)
  call blockDesc%deltas(del)
  dhdt = minval(del(1:NDIM))/dt

  !~ hy_fluxBuf[XYZ] represents (sum(F_fine) - F_coarse) on 
  !~ coarse side of f/c boundary, 0 elsewhere

  nullify(p_fluxBufX);nullify(p_fluxBufY);nullify(p_fluxBufZ)
  !These pointers allow us to use hy_fluxBuf[XYZ] (whose limits
  !are hard coded in FBS mode) within the varying loop limits below

  p_fluxBufX(1:NFLUXES,lo(IAXIS):hi(IAXIS)+1,lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) => hya_flxbx
#if NDIM>1  
  p_fluxBufY(1:NFLUXES,lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS)+1,lo(KAXIS):hi(KAXIS)) => hya_flxby
#if NDIM==3
  p_fluxBufZ(1:NFLUXES,lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)+1) => hya_flxbz
#endif
#endif

  if (hy_geometry /= CARTESIAN) then
     lev = blockDesc%level

     ! most of the following could and should be moved to geoFacs
     allocate(faceAreas(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
     call Grid_getCellFaceAreas(IAXIS,lev,loGC,hiGC,faceAreas)
     allocate(cellVolumes(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC))
     call Grid_getCellVolumes(lev,loGC,hiGC,cellVolumes)
     
     allocate(xCenter(xLoGC:xHiGC))
     allocate(xLeft(xLoGC:xHiGC))
     allocate(xRight(xLoGC:xHiGC)) 
     call Grid_getCellCoords(IAXIS, CENTER, lev, loGC, hiGC, xCenter)
     call Grid_getCellCoords(IAXIS, LEFT_EDGE, lev, loGC, hiGC, xLeft)
     call Grid_getCellCoords(IAXIS, RIGHT_EDGE, lev, loGC, hiGC, xRight)
  endif
 
!  !$omp parallel if (.FALSE.) &
!  !$omp default(none) &
!  !$omp firstprivate(vecLen)&
!  !$omp private(n,i,j,k,Vstar,facM,facP,range,eosData,massFraction,&
!  !$omp         emag,ekin,dx,dy,dz,fac,dFlux)&
!  !$omp shared(dt,solnData,hy_starState,p_fluxBufX,p_fluxBufY,p_fluxBufZ,hy_grav,&
!  !$omp        hya_flxbx,hya_flxby,hya_flxbz,xCenter,xLeft,xRight,blockDesc,&
!  !$omp        blkLimits,blkLimitsGC,hy_alphaGLM, hy_C_hyp,&
!  !$omp        dhdt, hy_smalldens, hy_smallE,del)

  ! Correct IAXIS sides
  !$omp do schedule(static) collapse(2)
  !Change limits to grownTile?
  do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
     do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        !! LOW side
        i = blkLimits(LOW,IAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        !~ (-) sign b/c fluxes are leaving on the low side
        dFlux = -p_fluxBufX(:,i  ,j  ,k  )
        dx = del(IAXIS)
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
        !call Eos_getData(range,vecLen,tempData,CENTER,eosData,massFraction)
        call Eos_getData(range,vecLen,solnData,CENTER,eosData,massFraction)
        !call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        !call Eos_putData(range,vecLen,tempData,CENTER,eosData)
        call Eos_putData(range,vecLen,solnData,CENTER,eosData)
        ! if (dFlux(HY_ENER) /= 0.) print *, 'a', Vstar(TEMP_VAR), dt/dx*fac*dFlux(HY_ENER)/(Vstar(ENER_VAR)*Vstar(DENS_VAR)), Vstar(VELX_VAR)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side
        i = blkLimits(HIGH,IAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        ! Positive b/c fluxes are entering on the high side
        dFlux = p_fluxBufX(:,i+1  ,j  ,k  )
        dx = del(IAXIS)
        ! Get geometric factors and sources
        call geoFacs(i,j,k,facM,facP)
        fac = facP
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i 
        range(LOW:HIGH,JAXIS) = j 
        range(LOW:HIGH,KAXIS) = k
       
        call Eos_getData(range,vecLen,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        call Eos_putData(range,vecLen,solnData,CENTER,eosData)

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
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -p_fluxBufY(:,i  ,j  ,k  )
        fac = 1.0
        dx = del(JAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i 
        range(LOW:HIGH,JAXIS) = j 
        range(LOW:HIGH,KAXIS) = k

        call Eos_getData(range,vecLen,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        call Eos_putData(range,vecLen,solnData,CENTER,eosData)
        ! Release pointers
        nullify(Vstar)

        !! HIGH side

        j = blkLimits(HIGH,JAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = p_fluxBufY(:,i  ,j+1  ,k  )
        fac = 1.0
        dx = del(JAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i 
        range(LOW:HIGH,JAXIS) = j 
        range(LOW:HIGH,KAXIS) = k
       
        call Eos_getData(range,vecLen,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        call Eos_putData(range,vecLen,solnData,CENTER,eosData)
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
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = -p_fluxBufZ(:,i  ,j  ,k  )
        fac = 1.0
        dx = del(KAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i
        range(LOW:HIGH,JAXIS) = j
        range(LOW:HIGH,KAXIS) = k
        
        call Eos_getData(range,vecLen,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        call Eos_putData(range,vecLen,solnData,CENTER,eosData)
        
        !release pointers
        nullify(Vstar)

        !! HIGH side
        k = blkLimits(HIGH,KAXIS)
        ! Point to old/intermediate states
        Vstar => solnData(:,i,j,k)
        ! Point to the correct fluxes
        dFlux = p_fluxBufZ(:,i  ,j  ,k+1  )
        fac = 1.0
        dx = del(KAXIS)
        call correctZone(Vstar,dFlux,dt,dx,fac)
        ! Update EOS
        !pos = (/i,j,k/)
        range(LOW:HIGH,IAXIS) = i
        range(LOW:HIGH,JAXIS) = j
        range(LOW:HIGH,KAXIS) = k
       
        call Eos_getData(range,vecLen,solnData,CENTER,eosData,massFraction)
        call Eos(MODE_DENS_EI,vecLen,eosData,massFraction)
        call Eos_putData(range,vecLen,solnData,CENTER,eosData) 
        ! Release pointers
        nullify(Vstar)
     enddo !j
  enddo !k
  !$omp end do
#endif
#endif

!  !$omp end parallel

  nullify(p_fluxBufX);nullify(p_fluxBufY);nullify(p_fluxBufZ)

  if (hy_geometry /= CARTESIAN) then
    deallocate(xCenter)
    deallocate(xLeft)
    deallocate(xRight)
    deallocate(faceAreas)
    deallocate(cellVolumes)
  endif

  deallocate(massFraction)
  deallocate(eosData)
  call blockDesc%releaseDataPtr(solnData,CENTER)
contains

  !!Update cell state based on flux deltas
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

    ! Now correct conserved vector with flux deltas
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
    facM = faceAreas(i  ,j,k)*dx/cellVolumes(i,j,k)
    facP = faceAreas(i+1,j,k)*dx/cellVolumes(i,j,k)

    if (xCenter(i) < 0.0) then
       facM = 0.
       facP = 0.
    end if
  end subroutine geoFacs

end subroutine hy_rk_correctFluxes
