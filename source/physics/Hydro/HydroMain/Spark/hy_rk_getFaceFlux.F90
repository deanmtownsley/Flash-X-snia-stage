!!****if* source/physics/Hydro/HydroMain/M1/hy_rk_getFaceFlux
!!
!!  NAME
!!
!!  hy_rk_getFaceFlux
!!
!!  SYNOPSIS
!!
!!  call hy_rk_getFaceFlux ( type(Grid_tile_t) :: blockDesc ) 
!!
!!  DESCRIPTION
!!  Initially stores grid data by 'pencils' (see subroutine setPencil() below),
!!  calls outside reconstruction subroutine, calls outisde Riemann solver, and 
!!  ultimately yields face fluxes for relevant directions.  These fluxes are lastly 
!!  saved to the grid multifabs (AMReX) for later access.
!!
!!  ARGUMENTS
!!  blockDesc-block/tile descriptor
!!  limits-sets limits of the fluxes.  Note this is a modified form typically yielded
!!         by blockDesc%limits b/c stage 1 of RK integration accounts for NSTENCIL 
!!         number of guard cells beyond the block interior.  Stage 2 only updates the 
!!         interior cells.
!!
!!***
!!!*** Already written in form of small function calls, so no need to insert kernel delimiters
!!Reorder(4):hy_starState,hy_fl[xyz] 
!!$subroutine hy_rk_getFaceFlux (blockDesc, limits)
!!$
!!$  use Grid_tile, ONLY : Grid_tile_t
!!$  use Hydro_data, ONLY : hy_threadWithinBlock, &
!!$                         hy_starState, hy_grav, hy_flattening, hy_flx, hy_fly, hy_flz
!!$  use Timers_interface, ONLY : Timers_start, Timers_stop
!!$
!!$  implicit none
!!$
!!$#include "Simulation.h"
!!$#include "constants.h"
!!$#include "Spark.h"
!!$#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS
!!$
!!$  type(Grid_tile_t)   :: blockDesc
!!$  integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
!!$
!!$!!!*** keyword here to indicate local allocation on device  
!!$
!!$  integer, dimension(LOW:HIGH,MDIM) :: blkLimits, blkLimitsGC
!!$
!!$  integer :: i1,i2,i3, n, g, dir, ierr,i
!!$
!!$  real, dimension(NFLUXES) :: flux
!!$
!!$  !integer, dimension(MDIM) :: datasize !deprecated
!!$  real, dimension(MDIM) :: del
!!$
!!$  ! 1D pencil of data
!!$  integer :: size1d !deprecated?
!!$  real, allocatable :: pencil(:,:)
!!$  integer, dimension(LOW:HIGH,MDIM) :: dirLims
!!$  real, dimension(NRECON) :: leftState, rightState
!!$  real, dimension(NRECON) :: uPlus, uMinus
!!$  real :: speed
!!$  
!!$!!!*** keyword here to indicate local allocation on device  
!!$  real, allocatable, target :: grv(:), shck(:), flat(:),flat3d(:,:,:)
!!$  real, pointer :: pgrv(:),pshck(:),pflat(:)!?,pflat3d(:,:,:)
!!$  logical :: inShock
!!$  integer :: pLo,pHi !low and high indices for pencil arrays
!!$  integer :: xLo,yLo,zLo,xHi,yHi,zHi,xLoGC,yLoGC,zLoGC,xHiGC,yHiGC,zHiGC
!!$
!!$  blkLimits(:,:)   = blockDesc%limits
!!$  blkLimitsGC(:,:) = blockDesc%blkLimitsGC
!!$
!!$  !convenience indices
!!$  xLo = blkLimits(LOW,IAXIS); xHi = blkLimits(HIGH,IAXIS)
!!$  yLo = blkLimits(LOW,JAXIS); yHi = blkLimits(HIGH,JAXIS)
!!$  zLo = blkLimits(LOW,KAXIS); zHi = blkLimits(HIGH,KAXIS)
!!$
!!$  xLoGC = blkLimitsGC(LOW,IAXIS); xHiGC = blkLimitsGC(HIGH,IAXIS)
!!$  yLoGC = blkLimitsGC(LOW,JAXIS); yHiGC = blkLimitsGC(HIGH,JAXIS)
!!$  zLoGC = blkLimitsGC(LOW,KAXIS); zHiGC = blkLimitsGC(HIGH,KAXIS)
!!$
!!$!!!*** Is this going to be locally allocated on the device?   
!!$  allocate(flat3d(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)) 
!!$
!!$  !set pencil indices based on dimension with largest spread
!!$  pLo = blkLimitsGC(LOW,maxloc(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM),1))
!!$  pHi = blkLimitsGC(HIGH,maxloc(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM),1))
!!$  
!!$  call blockDesc%deltas(del)
!!$  !datasize(1:MDIM) = blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM)+1
!!$  !size1d = maxval(datasize)
!!$  !size1d<-->pLo:pHi
!!$
!!$!!!*** Is this going to be locally allocated on the device?   
!!$!!!*** should all these allocations be taken out and handed to composer?
!!$  allocate(pencil(NRECON,pLo:pHi))
!!$  allocate(grv(pLo:pHi))
!!$  allocate(shck(pLo:pHi))
!!$  allocate(flat(pLo:pHi))
!!$  !allocate(temp(pLo:pHi))
!!$
!!$
!!$  if (hy_flattening) then
!!$!!!*** only needs transpiler
!!$     call flattening(flat3d, limits)
!!$  else
!!$     flat3d = 1.0
!!$  end if
!!$
!!$ !$omp parallel if(hy_threadWithinBlock .AND. NDIM > 1) &
!!$ !$omp default(none) &
!!$ !$omp private(i1,i2,i3,n,grv,shck,flat,inShock,flux,ierr,speed,&
!!$ !$omp         leftState,rightState,uPlus,uMinus,pencil,dir,dirLims,&
!!$ !$omp         pshck, pgrv,pflat)&
!!$ !$omp shared(hy_starState,blockDesc,pLo,pHi,xLoGC,yLoGC,zLoGC,&
!!$ !$omp        xHiGC,yHiGC,zHiGC,xLo,yLo,zLo,xHi,yHi,zHi,del,flat3d,limits)
!!$  !  Begin loop over zones
!!$  do dir = 1, NDIM
!!$     call setLoop(dir, dirLims)
!!$     
!!$     !$omp do schedule(guided) collapse(2)
!!$     do i3 = dirLims(LOW,3), dirLims(HIGH,3)
!!$        do i2 = dirLims(LOW,2), dirLims(HIGH,2)
!!$           !Define appropriate changing indices
!!$           select case(dir)
!!$           case (IAXIS)
!!$              pLo = xLoGC;pHi = xHiGC
!!$           case (JAXIS)
!!$              pLo = yLoGC;pHi = yHiGC
!!$           case (KAXIS)
!!$              pLo = zLoGC;pHi = zHiGC
!!$           end select
!!$
!!$           !pointers for purpose of reindexing
!!$           pshck(pLo:pHi) => shck
!!$           pgrv(pLo:pHi) => grv
!!$           pflat(pLo:pHi) => flat
!!$           
!!$!!!*** should be offloaded to composer or transpiler for sure
!!$           call setPencil(pencil,pgrv,pshck,pflat,i2,i3,dir)
!!$           !call setHSE(pencil,grv,dirLims(LOW,1)-1,del(dir))
!!$           i1 = dirLims(LOW,1) - 1 
!!$!!!*** needs only transpiler for the default implementation, still to look at the other ones
!!$           call reconstruct(uPlus, uMinus, &
!!$                pencil, pflat(i1), pLo, pHi, i1, del(dir))
!!$           !call resetHSE(pencil,grv,dirLims(LOW,1)-1,del(dir))
!!$           do i1 = dirLims(LOW,1), dirLims(HIGH,1)+1
!!$              ! cycle left-right state vectors
!!$              leftState = uPlus
!!$              ! Now do the reconstruction for ALL variables in one call
!!$              !call setHSE(pencil,grv,i1,del(dir))
!!$              call reconstruct(uPlus, uMinus, &
!!$                   pencil, pflat(i1), pLo, pHi, i1, del(dir))
!!$
!!$               !call resetHSE(pencil,grv,i1,del(dir))
!!$               rightState = uMinus
!!$               
!!$               ! Make sure that reconstruction has not introduced unphysical states
!!$               call ensurePhysicalState(leftState)
!!$               call ensurePhysicalState(rightState)
!!$               
!!$               ! Check for shocks in the zones involved in flux calculation
!!$               inShock = any(pshck(i1-1:i1) /= 0.0)
!!$               ! Now call the Riemann solver to compute fluxes
!!$               !!!*** maybe refactor a little bit, if primitive to conservative conversion should be done per block etc
!!$               !!!*** quite likely it is fine with just being inlined by the transpiler              
!!$               call riemann(dir,leftState(1:HY_NUM_VARS),&
!!$               rightState(1:HY_NUM_VARS),inShock,&
!!$               flux(1:HY_NUM_FLUX),speed,ierr)
!!$              ! Add artificial viscosity for strong-shock capturing
!!$!!!*** only transpiler needed
!!$              call avisc(flux,pencil,i1,NRECON,dir)
!!$
!!$              ! Here, we compute the species and mass scalar
!!$              ! fluxes based on the density flux and the reconstructed
!!$              ! mass scalar interface values
!!$
!!$!!!*** only transpiler needed              
!!$              call mscalarFluxes(flux,leftState(HY_NUM_VARS+1:NRECON),&
!!$                   rightState(HY_NUM_VARS+1:NRECON))
!!$
!!$              ! ***************
!!$              ! Fluxes computed for one face of this zone
!!$              ! Save the fluxes
!!$              ! ***************
!!$              call saveFluxes(dir,i1,i2,i3,flux)
!!$           end do  ! i
!!$        !release pointers
!!$        nullify(pshck)
!!$        nullify(pgrv)
!!$        nullify(pflat)
!!$        end do ! j
!!$     end do ! k
!!$     !$omp end do
!!$  end do ! dir
!!$ !$omp end parallel
!!$
!!$  deallocate(pencil)
!!$  deallocate(grv)
!!$  deallocate(shck)
!!$  deallocate(flat)
!!$  deallocate(flat3d)
!!$
!!$contains
!!$
!!$  !!Set loop dimensions based on the direction of the pencil set
!!$  subroutine setLoop(dir, dirLims)
!!$    implicit none
!!$    integer, intent(IN) :: dir
!!$    integer, intent(OUT), dimension(LOW:HIGH,MDIM) :: dirLims
!!$    select case(dir)
!!$    case (IAXIS)
!!$       dirLims(:,:) = limits(:,:)
!!$    case (JAXIS)
!!$       dirLims(:,1) = limits(:,JAXIS)
!!$       dirLims(:,2) = limits(:,IAXIS)
!!$       dirLims(:,3) = limits(:,KAXIS)
!!$    case (KAXIS)
!!$       dirLims(:,1) = limits(:,KAXIS)
!!$       dirLims(:,2) = limits(:,IAXIS)
!!$       dirLims(:,3) = limits(:,JAXIS)
!!$    end select
!!$
!!$  end subroutine setLoop
!!$
!!$  !! the 'pencil' holds a 1D array of the solution data to be operated
!!$  !! on.  It is unrolled this way so that all data that are needed for
!!$  !! interpolation and flux calculation are truly contiguous in memory
!!$  !! space.  And so that it all fits in cache at once...
!!$  !! The maximum amount of calculation is done on these data prior to
!!$  !! reseting the pencil data to a new ray through the block.
!!$  !! If I were talented at ASCII art, I would make a diagram...
!!$  subroutine setPencil(pencil,grv,shck,flat,i2,i3,dir)
!!$    implicit none
!!$    real, intent(INOUT) :: pencil(HY_NUM_VARS+NSPECIES+NMASS_SCALARS,pLo:pHi)
!!$    real, intent(INOUT) :: grv(pLo:pHi)!size1d
!!$    real, intent(INOUT) :: shck(pLo:pHi)
!!$    real, intent(INOUT) :: flat(pLo:pHi)
!!$    integer, intent(IN) :: i2,i3,dir
!!$    integer :: n
!!$
!!$    select case(dir)
!!$    case (IAXIS)
!!$       pencil(HY_DENS,:) = hy_starState(DENS_VAR,:,i2,i3)
!!$       pencil(HY_VELX,:) = hy_starState(VELX_VAR,:,i2,i3)
!!$       pencil(HY_VELY,:) = hy_starState(VELY_VAR,:,i2,i3)
!!$       pencil(HY_VELZ,:) = hy_starState(VELZ_VAR,:,i2,i3)
!!$       pencil(HY_PRES,:) = hy_starState(PRES_VAR,:,i2,i3)
!!$       pencil(HY_GAMC,:) = hy_starState(GAMC_VAR,:,i2,i3)
!!$       pencil(HY_RHOE,:) = hy_starState(DENS_VAR,:,i2,i3)*hy_starState(EINT_VAR,:,i2,i3)
!!$       !!!*** should we consider replacing these ifdefs by keywords for the composer?
!!$#ifdef SPARK_GLM
!!$       pencil(HY_MAGX,:) = hy_starState(MAGX_VAR,:,i2,i3)
!!$       pencil(HY_MAGY,:) = hy_starState(MAGY_VAR,:,i2,i3)
!!$       pencil(HY_MAGZ,:) = hy_starState(MAGZ_VAR,:,i2,i3)
!!$       pencil(HY_PSIB,:) = hy_starState(PSIB_VAR,:,i2,i3)
!!$#endif
!!$#if NSPECIES+NMASS_SCALARS>0
!!$       do n=SPECIES_BEGIN, MASS_SCALARS_END
!!$          pencil(HY_NUM_VARS+1+n-SPECIES_BEGIN,:)    = hy_starState(n,:,i2,i3)
!!$       enddo
!!$#endif
!!$#ifdef GRAVITY
!!$#ifdef GPOT_VAR
!!$       grv(:) = hy_starState(GPOT_VAR,:,i2,i3)
!!$#else
!!$       grv(:) = hy_grav(IAXIS,:,i2,i3)
!!$#endif
!!$#endif
!!$#ifdef SHOK_VAR
!!$       shck(:) = hy_starState(SHOK_VAR,:,i2,i3)
!!$#else
!!$       shck(:) = 0.0
!!$#endif
!!$       flat(:) = flat3d(:,i2,i3)
!!$    case (JAXIS)
!!$       pencil(HY_DENS,:) = hy_starState(DENS_VAR,i2,:,i3)
!!$       pencil(HY_VELX,:) = hy_starState(VELX_VAR,i2,:,i3)
!!$       pencil(HY_VELY,:) = hy_starState(VELY_VAR,i2,:,i3)
!!$       pencil(HY_VELZ,:) = hy_starState(VELZ_VAR,i2,:,i3)
!!$       pencil(HY_PRES,:) = hy_starState(PRES_VAR,i2,:,i3)
!!$       pencil(HY_GAMC,:) = hy_starState(GAMC_VAR,i2,:,i3)
!!$       pencil(HY_RHOE,:) = hy_starState(DENS_VAR,i2,:,i3)*hy_starState(EINT_VAR,i2,:,i3)
!!$#ifdef SPARK_GLM
!!$       pencil(HY_MAGX,:) = hy_starState(MAGX_VAR,i2,:,i3)
!!$       pencil(HY_MAGY,:) = hy_starState(MAGY_VAR,i2,:,i3)
!!$       pencil(HY_MAGZ,:) = hy_starState(MAGZ_VAR,i2,:,i3)
!!$       pencil(HY_PSIB,:) = hy_starState(PSIB_VAR,i2,:,i3)
!!$#endif
!!$#if NSPECIES+NMASS_SCALARS>0
!!$       do n=SPECIES_BEGIN, MASS_SCALARS_END
!!$          pencil(HY_NUM_VARS+1+n-SPECIES_BEGIN,:)    = hy_starState(n,i2,:,i3)
!!$       enddo
!!$#endif
!!$#ifdef GRAVITY
!!$#ifdef GPOT_VAR
!!$       grv(:) = hy_starState(GPOT_VAR,i2,:,i3)
!!$#else
!!$       grv(:) = hy_grav(JAXIS,i2,:,i3)
!!$#endif
!!$#endif
!!$#ifdef SHOK_VAR
!!$       shck(:) = hy_starState(SHOK_VAR,i2,:,i3)
!!$#else
!!$       shck(:) = 0.0
!!$#endif
!!$       flat(:) = flat3d(i2,:,i3)
!!$    case (KAXIS)
!!$       pencil(HY_DENS,:) = hy_starState(DENS_VAR,i2,i3,:)
!!$       pencil(HY_VELX,:) = hy_starState(VELX_VAR,i2,i3,:)
!!$       pencil(HY_VELY,:) = hy_starState(VELY_VAR,i2,i3,:)
!!$       pencil(HY_VELZ,:) = hy_starState(VELZ_VAR,i2,i3,:)
!!$       pencil(HY_PRES,:) = hy_starState(PRES_VAR,i2,i3,:)
!!$       pencil(HY_GAMC,:) = hy_starState(GAMC_VAR,i2,i3,:)
!!$       pencil(HY_RHOE,:) = hy_starState(DENS_VAR,i2,i3,:)*hy_starState(EINT_VAR,i2,i3,:)
!!$#ifdef SPARK_GLM
!!$       pencil(HY_MAGX,:) = hy_starState(MAGX_VAR,i2,i3,:)
!!$       pencil(HY_MAGY,:) = hy_starState(MAGY_VAR,i2,i3,:)
!!$       pencil(HY_MAGZ,:) = hy_starState(MAGZ_VAR,i2,i3,:)
!!$       pencil(HY_PSIB,:) = hy_starState(PSIB_VAR,i2,i3,:)
!!$#endif
!!$#if NSPECIES+NMASS_SCALARS>0
!!$       do n=SPECIES_BEGIN, MASS_SCALARS_END
!!$          pencil(HY_NUM_VARS+1+n-SPECIES_BEGIN,:)    = hy_starState(n,i2,i3,:)
!!$       enddo
!!$#endif
!!$#ifdef GRAVITY
!!$#ifdef GPOT_VAR
!!$       grv(:) = hy_starState(GPOT_VAR,i2,i3,:)
!!$#else
!!$       grv(:) = hy_grav(KAXIS,i2,i3,:)
!!$#endif
!!$#endif
!!$#ifdef SHOK_VAR
!!$       shck(:) = hy_starState(SHOK_VAR,i2,i3,:)
!!$#else
!!$       shck(:) = 0.0
!!$#endif
!!$       flat(:) = flat3d(i2,i3,:)
!!$    end select
!!$  end subroutine setPencil
!!$
!!$
!!$  subroutine saveFluxes(dir,i1,i2,i3,flux)
!!$    integer, intent(IN) :: dir,i1,i2,i3
!!$    real, dimension(NFLUXES), intent(IN) :: flux
!!$    select case(dir)
!!$    case(IAXIS)
!!$     hy_flx(:,i1,i2,i3) = flux(:)
!!$
!!$    case (JAXIS)
!!$     hy_fly(:,i2,i1,i3) = flux(:)
!!$
!!$    case (KAXIS)
!!$     hy_flz(:,i2,i3,i1) = flux(:)
!!$
!!$    end select
!!$  end subroutine saveFluxes
!!$
!!$  !!Account for fluxes that are proportional to mass flux 
!!$  subroutine mscalarFluxes(flux,XL,XR)
!!$    implicit none
!!$    real, intent(INOUT) :: flux(NFLUXES)
!!$    real, target, dimension(NSPECIES+NMASS_SCALARS), intent(IN) :: XL, XR
!!$    real, pointer :: Xstar(:)
!!$#if NSPECIES+NMASS_SCALARS==0
!!$    return
!!$#else
!!$    if (flux(HY_MASS) > 0.) then
!!$       Xstar => XL
!!$    else
!!$       Xstar => XR
!!$    endif
!!$    flux(HY_NUM_FLUX+1:NFLUXES) = Xstar*flux(HY_MASS)
!!$#endif
!!$  end subroutine mscalarFluxes
!!$
!!$  subroutine setHSE(V,grv,i1,dx)
!!$    implicit none
!!$    real, intent(INOUT) :: V(HY_NUM_VARS+NSPECIES+NMASS_SCALARS,pLo:pHi)!size1d
!!$    real, intent(IN) :: grv(pLo:pHi)!size1d
!!$    integer, intent(IN) :: i1
!!$    real, intent(IN) :: dx
!!$    real :: accelM, accelP
!!$#ifndef GRAVITY
!!$    return
!!$#else
!!$#ifdef GPOT_VAR
!!$    accelM = grv(i1)   - grv(i1-1)
!!$    accelP = grv(i1+1) - grv(i1)
!!$#else
!!$    accelM = 0.5*(grv(i1-1) + grv(i1)  )*dx
!!$    accelP = 0.5*(grv(i1)   + grv(i1+1))*dx
!!$#endif
!!$    V(HY_PRES,i1-1) = V(HY_PRES,i1-1) - &
!!$         0.5*(V(HY_DENS,i1-1)+V(HY_DENS,i1)) * accelM
!!$    V(HY_PRES,i1+1) = V(HY_PRES,i1+1) + &
!!$         0.5*(V(HY_DENS,i1)+V(HY_DENS,i1+1)) * accelP
!!$#endif /* GRAVITY */
!!$  end subroutine setHSE
!!$
!!$  subroutine resetHSE(V,grv,i1,dx)
!!$    implicit none
!!$    real, intent(INOUT) :: V(HY_NUM_VARS+NSPECIES+NMASS_SCALARS,pLo:pHi)
!!$    real, intent(IN) :: grv(pLo:pHi)!size1d
!!$    integer, intent(IN) :: i1
!!$    real, intent(IN) :: dx
!!$    real :: accelM, accelP
!!$#ifndef GRAVITY
!!$    return
!!$#else
!!$#ifdef GPOT_VAR
!!$    accelM = grv(i1)   - grv(i1-1)
!!$    accelP = grv(i1+1) - grv(i1)
!!$#else
!!$    accelM = 0.5*(grv(i1-1) + grv(i1)  )*dx
!!$    accelP = 0.5*(grv(i1)   + grv(i1+1))*dx
!!$#endif
!!$    V(HY_PRES,i1-1) = V(HY_PRES,i1-1) + &
!!$         0.5*(V(HY_DENS,i1-1)+V(HY_DENS,i1)) * accelM
!!$    V(HY_PRES,i1+1) = V(HY_PRES,i1+1) - &
!!$         0.5*(V(HY_DENS,i1)+V(HY_DENS,i1+1)) * accelP
!!$#endif /* GRAVITY */
!!$  end subroutine resetHSE
!!$
!!$  !! Artificial viscosity as in Colella and Woodward
!!$  !! Calculation of cvisc here is different from the UHD solver
!!$  !! I am not including the transverse velocity terms to keep everything
!!$  !! neat and 1D for the pencil data
!!$  subroutine avisc(flux,V,i1,nvars,dir)
!!$    use Hydro_data, ONLY : hy_cvisc
!!$    implicit none
!!$    real, intent(INOUT) :: flux(NFLUXES)
!!$    integer, intent(IN) :: nvars
!!$    real, intent(IN) :: V(nvars,pLo:pHi)!size1d
!!$    integer, intent(IN) :: i1,dir
!!$    real :: cvisc, VenerLo, VenerHi
!!$
!!$    cvisc = hy_cvisc*max(-(V(HY_VELX+dir-1,i1) - V(HY_VELX+dir-1,i1-1)),0.)
!!$
!!$    ! Construct minus and plus TOTAL energy densities
!!$    VenerLo = V(HY_DENS,i1-1)*0.5*(dot_product(V(HY_VELX:HY_VELZ,i1-1),V(HY_VELX:HY_VELZ,i1-1)))&
!!$            + V(HY_RHOE,i1-1)
!!$    VenerHi = V(HY_DENS,i1)*0.5*(dot_product(V(HY_VELX:HY_VELZ,i1),V(HY_VELX:HY_VELZ,i1)))&
!!$            + V(HY_RHOE,i1)
!!$
!!$    flux(HY_MASS:HY_ENER) = &
!!$         flux(HY_MASS:HY_ENER) &
!!$         +cvisc*(/V(HY_DENS,i1-1)                 - V(HY_DENS,i1)&
!!$         ,        V(HY_DENS,i1-1)*V(HY_VELX,i1-1) - V(HY_DENS,i1)*V(HY_VELX,i1)&
!!$         ,        V(HY_DENS,i1-1)*V(HY_VELY,i1-1) - V(HY_DENS,i1)*V(HY_VELY,i1)&
!!$         ,        V(HY_DENS,i1-1)*V(HY_VELZ,i1-1) - V(HY_DENS,i1)*V(HY_VELZ,i1)&
!!$         ,        VenerLo                         - VenerHi/)
!!$#ifdef SPARK_GLM
!!$    flux(HY_FMGX:HY_FPSI) = &
!!$         flux(HY_FMGX:HY_FPSI) &
!!$         +cvisc*(/V(HY_MAGX,i1-1)                 - V(HY_MAGX,i1)&
!!$         ,        V(HY_MAGY,i1-1)                 - V(HY_MAGY,i1)&
!!$         ,        V(HY_MAGZ,i1-1)                 - V(HY_MAGZ,i1)&
!!$         ,        V(HY_PSIB,i1-1)                 - V(HY_PSIB,i1)/)
!!$#endif
!!$  end subroutine avisc
!!$
!!$  !~ Flattening has not been tested yet in FLASH5, only 1D & 2D runs so far.
!!$  subroutine flattening(flat3d,limits)
!!$    !! This follows Miller & Colella 2002
!!$    use Hydro_data, ONLY : hy_starState
!!$    implicit none
!!$    integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
!!$    !real, intent(OUT) :: flat3d(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
!!$    !real :: flatTilde(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
!!$    real, intent(OUT) :: flat3d(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)
!!$    real :: flatTilde(NDIM,xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)
!!$    real :: beta, Z
!!$    real, parameter :: betaMin = 0.75, betaMax = 0.85
!!$    real, parameter :: Zmin = 0.25, Zmax = 0.75
!!$    integer :: i,j,k, kx, ky, kz
!!$    kx = 1
!!$    ky = 0
!!$    kz = 0
!!$#if NDIM>1
!!$    ky = 1
!!$#if NDIM==3
!!$    kz = 1
!!$#endif
!!$#endif
!!$    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
!!$       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
!!$          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
!!$             !1e-99 used to be TINY(1.0) but that produced Seg Faults
!!$             beta = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
!!$                  / max(1e-99, abs(hy_starState(PRES_VAR,i+2,j,k)-hy_starState(PRES_VAR,i-2,j,k)))
!!$             Z    = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
!!$                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!!$             flatTilde(IAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!!$             if (hy_starState(VELX_VAR,i+1,j,k)<hy_starState(VELX_VAR,i,j,k)) then
!!$                flatTilde(IAXIS,i,j,k) = max(flatTilde(IAXIS,i,j,k), &
!!$                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!!$             else
!!$                flatTilde(IAXIS,i,j,k) = 1.0
!!$             end if
!!$#if NDIM>1
!!$             beta = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
!!$                  / max(1.e-99,abs(hy_starState(PRES_VAR,i,j+2,k)-hy_starState(PRES_VAR,i,j-2,k)))
!!$             Z    = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
!!$                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!!$             flatTilde(JAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!!$             if (hy_starState(VELY_VAR,i,j+1,k)<hy_starState(VELY_VAR,i,j,k)) then
!!$                flatTilde(JAXIS,i,j,k) = max(flatTilde(JAXIS,i,j,k), &
!!$                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!!$             else
!!$                flatTilde(JAXIS,i,j,k) = 1.0
!!$             end if
!!$#if NDIM==3
!!$             beta = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
!!$                  / max(1.e-99,abs(hy_starState(PRES_VAR,i,j,k+2)-hy_starState(PRES_VAR,i,j,k-2)))
!!$             Z    = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
!!$                  / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!!$             flatTilde(KAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!!$             if (hy_starState(VELZ_VAR,i,j,k+1)<hy_starState(VELZ_VAR,i,j,k)) then
!!$                flatTilde(KAXIS,i,j,k) = max(flatTilde(KAXIS,i,j,k), &
!!$                     &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!!$             else
!!$                flatTilde(KAXIS,i,j,k) = 1.0
!!$             end if
!!$#endif
!!$#endif
!!$          end do
!!$       end do
!!$    end do
!!$    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
!!$       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
!!$          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
!!$             flat3d(i,j,k) = minval(flatTilde(1:NDIM,i,j,k))
!!$#ifdef FLAT_VAR
!!$             hy_starState(FLAT_VAR,i,j,k) = flat3d(i,j,k)
!!$#endif
!!$          end do
!!$       end do
!!$    end do
!!$  end subroutine flattening
!!$
!!$  subroutine ensurePhysicalState(state)
!!$    use Hydro_data, ONLY : hy_smalldens, hy_smallE, hy_smallpres, hy_smallX
!!$    implicit none
!!$    real, dimension(NRECON), intent(INOUT), target  :: state
!!$    integer :: s
!!$    real :: spcSumInv
!!$    real, pointer :: spc(:)
!!$    state(HY_DENS) = max(hy_smalldens, state(HY_DENS))
!!$    state(HY_PRES) = max(hy_smallpres, state(HY_PRES))
!!$#if NSPECIES>0
!!$    ! Limit and renormalize the species.
!!$    spc => state(HY_NUM_VARS+1:HY_NUM_VARS+NSPECIES)
!!$    do s = 1, NSPECIES
!!$       spc(s) = max(hy_smallX,min(1.0,spc(s)))
!!$    end do
!!$    spcSumInv = 1./sum(spc(1:NSPECIES))
!!$    spc = spc*spcSumInv
!!$#endif
!!$  end subroutine ensurePhysicalState
!!$
!!$
!!$end subroutine hy_rk_getFaceFlux



!!****if* source/physics/Hydro/HydroMain/M1/hy_rk_getFaceFlux
!!
!!  NAME
!!
!!  hy_rk_getFaceFlux_offloaded
!!
!!  SYNOPSIS
!!
!!  call hy_rk_getFaceFlux_offloaded ( type(Grid_tile_t) :: blockDesc )
!!
!!  DESCRIPTION
!!  Initially stores grid data by 'pencils' (see subroutine setPencil() below),
!!  calls outside reconstruction subroutine, calls outisde Riemann solver, and 
!!  ultimately yields face fluxes for relevant directions.  These fluxes are lastly 
!!  saved to the grid multifabs (AMReX) for later access.
!!
!!  ARGUMENTS
!!  blockDesc-block/tile descriptor
!!  limits-sets limits of the fluxes.  Note this is a modified form typically yielded
!!         by blockDesc%limits b/c stage 1 of RK integration accounts for NSTENCIL 
!!         number of guard cells beyond the block interior.  Stage 2 only updates the 
!!         interior cells.
!!
!!***
!!!*** Already written in form of small function calls, so no need to insert kernel delimiters
!!Reorder(4):hy_starState,hy_fl[xyz]
subroutine hy_rk_getFaceFlux(limits,blkLimits,blkLimitsGC)

   use Hydro_data, ONLY : hy_threadWithinBlock, &
                          hy_starState, hy_grav, hy_flattening, hy_flx, hy_fly, hy_flz, &
                          snake, uPlusArray, uMinusArray, flat, grv, shck, flux, &
                          flat3d
   use Timers_interface, ONLY : Timers_start, Timers_stop
   use Driver_interface, ONLY : Driver_abort
 
   implicit none
 
#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS
 
   integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits,blkLimits, blkLimitsGC
 
   integer :: i1,i2,i3, n, g, dir, ierr
   integer, dimension(3) :: guardCells
   integer, dimension(LOW:HIGH,MDIM) :: dirLims
   integer :: i_global, j_global, k_global,v
   character(len = 2) :: dir_str

    
   !$omp target data map(to: dir, dirLims, guardCells)

   if (hy_flattening) then
      !!!*** only needs transpiler
      call Driver_abort( "Flattening has not be implemented with GPU offloading nor has it been tested in FLASH5.")
      ! call flattening(flat3d, limits)
   else
      ! call Timers_start("flat3d")
      !$omp target teams distribute parallel do collapse(3) shared(blkLimitsGC,flat3d) private(i1,i2,i3) default(none) map(to:blkLimitsGC)!! TODO: Set this once for both rk steps.
      do i3 = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
         do i2 = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
            do i1 = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
               flat3d(i1,i2,i3) = 1.0
            enddo
         enddo
      enddo
      ! call Timers_stop("flat3d")
   end if



   ! !  Begin loop over zones
   do dir = 1,  NDIM

      select case(dir)
      case (IAXIS)
         dirLims(:,:) = limits(:,:)
         guardCells(1) = (limits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
         guardCells(2) = (limits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
         guardCells(3) = (limits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))
         dir_str = "_x"
      case (JAXIS)
         dirLims(:,1) = limits(:,JAXIS)
         dirLims(:,2) = limits(:,IAXIS)
         dirLims(:,3) = limits(:,KAXIS)
         guardCells(1) = (limits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
         guardCells(2) = (limits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
         guardCells(3) = (limits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))
         dir_str = "_y"
      case (KAXIS)
         dirLims(:,1) = limits(:,KAXIS)
         dirLims(:,2) = limits(:,IAXIS)
         dirLims(:,3) = limits(:,JAXIS)
         guardCells(1) = (limits(LOW,KAXIS) - blkLimitsGC(LOW,KAXIS))
         guardCells(2) = (limits(LOW,IAXIS) - blkLimitsGC(LOW,IAXIS))
         guardCells(3) = (limits(LOW,JAXIS) - blkLimitsGC(LOW,JAXIS))
         dir_str = "_z"
      end select 
      !$omp target update to(dir,guardCells,dirLims)
      ! call Timers_start("snake_init"//dir_str) 

      call initializeSnake(dir, dirLims, guardCells) ! This is offloaded internally 


     !$omp target teams distribute parallel do collapse(4) &
     !$omp  private(i1,i2,i3,v) shared(dirLims,guardCells,uPlusArray, uMinusArray, snake, flat) default(none)
     do i3 =1 + guardCells(3) , 1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
        do i2 =1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
           do i1 =1 + guardCells(1) - 1 , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
              do v=1,NRECON
                  call reconstruct(i1,i2,i3,v)
              enddo 
           enddo
        enddo
     enddo
   !   call reconstructOffloaded_split_inlined(dirLims,guardCells)
      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3,v) shared(dirLims,dir,guardCells)
      do i3 = 1 + guardCells(3) , 1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) - 1 , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call ensurePhysicalState(uPlusArray(:,i1,i2,i3))
               call ensurePhysicalState(uMinusArray(:,i1,i2,i3))
            enddo 
         enddo
      enddo


      ! call Timers_start("riemann"//dir_str)
      !  Now call the Riemann solver to compute fluxes
      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells,flux) default(none)
       do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
          do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
             do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
                 call riemann(i1,i2,i3,dir)   
            enddo
         enddo
      enddo
      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells,flux) default(none)
      do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call avisc(i1,i2,i3,dir,NRECON)
               call mscalarFluxes(i1,i2,i3)
           enddo
        enddo
     enddo

      ! call Timers_stop("riemann"//dir_str)

       ! ***************
       ! Fluxes computed for one face of this zone
       ! Save the fluxes
       ! ***************
      !  call Timers_start("saveFluxes"//dir_str)

      !$omp target teams distribute parallel do collapse(3) &
      !$omp private(i1,i2,i3) shared(dirLims,dir,guardCells) default(none)
       do i3 = 1 + guardCells(3) ,1 + guardCells(3) + dirLims(HIGH,3) - dirLims(LOW,3)
         do i2 = 1 + guardCells(2) , 1 + guardCells(2) + dirLims(HIGH,2) - dirLims(LOW,2)
            do i1 = 1 + guardCells(1) , 1 + guardCells(1) + dirLims(HIGH,1) - dirLims(LOW,1) + 1
               call saveFluxes(dir,i1,i2,i3,guardCells,dirLims)
             enddo
          enddo
       enddo 
      !  call Timers_stop("saveFluxes"//dir_str)
      
    enddo
   !$omp end target data 
 contains
 
   !!Set loop dimensions based on the direction of the pencil set
   subroutine setLoop(dir)
      use Hydro_data, ONLY: dirLims
     implicit none
     integer, intent(IN) :: dir
     !$omp declare target
     select case(dir)
     case (IAXIS)
        dirLims(:,:) = limits(:,:)
     case (JAXIS)
        dirLims(:,1) = limits(:,JAXIS)
        dirLims(:,2) = limits(:,IAXIS)
        dirLims(:,3) = limits(:,KAXIS)
     case (KAXIS)
        dirLims(:,1) = limits(:,KAXIS)
        dirLims(:,2) = limits(:,IAXIS)
        dirLims(:,3) = limits(:,JAXIS)
     end select
 
   end subroutine setLoop
 
   subroutine saveFluxes(dir,i1,i2,i3,guardCells,dirLims)
      USE Hydro_data, ONLY: flux, hy_flx, hy_fly, hy_flz
     integer, intent(IN) :: dir,i1,i2,i3
     integer, dimension(3) :: guardCells
     integer, dimension(LOW:HIGH,MDIM) :: dirLims
     integer :: i_global, j_global, k_global
     !$omp declare target
     i_global = -1 + i1 + dirLims(LOW,1) - guardCells(1)
     j_global = -1 + i2 + dirLims(LOW,2) - guardCells(2)
     k_global = -1 + i3 + dirLims(LOW,3) - guardCells(3)
     select case(dir)
     case(IAXIS)
        hy_flx(:,i_global,j_global,k_global) = flux(:,i1,i2,i3)
     case (JAXIS)
        hy_fly(:,j_global,i_global,k_global) = flux(:,i1,i2,i3)
     case (KAXIS)
        hy_flz(:,j_global,k_global,i_global) = flux(:,i1,i2,i3)
     end select
   end subroutine saveFluxes
 
!!Account for fluxes that are proportional to mass flux 
subroutine mscalarFluxes(i1,i2,i3)
   USE Hydro_data, ONLY: flux, uPlusArray, uMinusArray
   implicit none
   real, pointer :: Xstar(:)
   integer, intent(IN) :: i1,i2,i3
   !$omp declare target
#if NSPECIES+NMASS_SCALARS==0
   return
#else
   if (flux(HY_MASS,i1,i2,i3) > 0.) then
      Xstar => uPlusArray(HY_NUM_VARS+1:NRECON,i1-1,i2,i3)
   else
      Xstar => uMinusArray(HY_NUM_VARS+1:NRECON,i1,i2,i3)
   endif
   flux(HY_NUM_FLUX+1:NFLUXES,i1,i2,i3) = Xstar*flux(HY_MASS,i1,i2,i3)
#endif
end subroutine mscalarFluxes
 
 
!! Artificial viscosity as in Colella and Woodward
!! Calculation of cvisc here is different from the UHD solver
!! I am not including the transverse velocity terms to keep everything
!! neat and 1D for the pencil data
subroutine avisc(i1,i2,i3,dir,nvars) !flux,V,i1,nvars,dir)
   use Hydro_data, ONLY : hy_cvisc, flux, snake
   implicit none
   integer, intent(IN) :: nvars
   integer, intent(IN) :: i1,i2,i3,dir
   real :: cvisc, VenerLo, VenerHi
   !$omp declare target
   cvisc = hy_cvisc*max(-(snake(HY_VELX+dir-1,i1,i2,i3) - snake(HY_VELX+dir-1,i1-1,i2,i3)),0.)

   ! Construct minus and plus TOTAL energy densities
   VenerLo = snake(HY_DENS,i1-1,i2,i3)*0.5*(dot_product(snake(HY_VELX:HY_VELZ,i1-1,i2,i3),snake(HY_VELX:HY_VELZ,i1-1,i2,i3)))&
            + snake(HY_RHOE,i1-1,i2,i3)
   VenerHi = snake(HY_DENS,i1,i2,i3)*0.5*(dot_product(snake(HY_VELX:HY_VELZ,i1,i2,i3),snake(HY_VELX:HY_VELZ,i1,i2,i3)))&
            + snake(HY_RHOE,i1,i2,i3)

   flux(HY_MASS:HY_ENER,i1,i2,i3) = &
         flux(HY_MASS:HY_ENER,i1,i2,i3) &
         +cvisc*(/snake(HY_DENS,i1-1,i2,i3)                 - snake(HY_DENS,i1,i2,i3)&
         ,        snake(HY_DENS,i1-1,i2,i3)*snake(HY_VELX,i1-1,i2,i3) - snake(HY_DENS,i1,i2,i3)*snake(HY_VELX,i1,i2,i3)&
         ,        snake(HY_DENS,i1-1,i2,i3)*snake(HY_VELY,i1-1,i2,i3) - snake(HY_DENS,i1,i2,i3)*snake(HY_VELY,i1,i2,i3)&
         ,        snake(HY_DENS,i1-1,i2,i3)*snake(HY_VELZ,i1-1,i2,i3) - snake(HY_DENS,i1,i2,i3)*snake(HY_VELZ,i1,i2,i3)&
         ,        VenerLo                         - VenerHi/)
#ifdef SPARK_GLM
   flux(HY_FMGX:HY_FPSI,i1,i2,i3) = &
         flux(HY_FMGX:HY_FPSI,i1,i2,i3) &
         +cvisc*(/snake(HY_MAGX,i1-1,i2,i3)                 - snake(HY_MAGX,i1,i2,i3)&
         ,        snake(HY_MAGY,i1-1,i2,i3)                 - snake(HY_MAGY,i1,i2,i3)&
         ,        snake(HY_MAGZ,i1-1,i2,i3)                 - snake(HY_MAGZ,i1,i2,i3)&
         ,        snake(HY_PSIB,i1-1,i2,i3)                 - snake(HY_PSIB,i1,i2,i3)/)
#endif
end subroutine avisc
 
! !~ Flattening has not been tested yet in FLASH5, only 1D & 2D runs so far.
! subroutine flattening(flat3d,limits)
!    !! This follows Miller & Colella 2002
!    use Hydro_data, ONLY : hy_starState
!    implicit none
!    integer, intent(IN), dimension(LOW:HIGH,MDIM) :: limits
!    !real, intent(OUT) :: flat3d(GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
!    !real :: flatTilde(NDIM,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
!    real, intent(OUT) :: flat3d(xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)
!    real :: flatTilde(NDIM,xLoGC:xHiGC,yLoGC:yHiGC,zLoGC:zHiGC)
!    real :: beta, Z
!    real, parameter :: betaMin = 0.75, betaMax = 0.85
!    real, parameter :: Zmin = 0.25, Zmax = 0.75
!    integer :: i,j,k, kx, ky, kz
!       pLo = blkLimitsGC(LOW,maxloc(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM),1))
!    pHi = blkLimitsGC(HIGH,maxloc(blkLimitsGC(HIGH,1:MDIM)-blkLimitsGC(LOW,1:MDIM),1))
!    kx = 1
!    ky = 0
!    kz = 0
! #if NDIM>1
!    ky = 1
! #if NDIM==3
!    kz = 1
! #endif
! #endif
!    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
!       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
!          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
!             !1e-99 used to be TINY(1.0) but that produced Seg Faults
!             beta = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
!                   / max(1e-99, abs(hy_starState(PRES_VAR,i+2,j,k)-hy_starState(PRES_VAR,i-2,j,k)))
!             Z    = abs(hy_starState(PRES_VAR,i+1,j,k)-hy_starState(PRES_VAR,i-1,j,k)) &
!                   / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!             flatTilde(IAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!             if (hy_starState(VELX_VAR,i+1,j,k)<hy_starState(VELX_VAR,i,j,k)) then
!                flatTilde(IAXIS,i,j,k) = max(flatTilde(IAXIS,i,j,k), &
!                      &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!             else
!                flatTilde(IAXIS,i,j,k) = 1.0
!             end if
! #if NDIM>1
!             beta = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
!                   / max(1.e-99,abs(hy_starState(PRES_VAR,i,j+2,k)-hy_starState(PRES_VAR,i,j-2,k)))
!             Z    = abs(hy_starState(PRES_VAR,i,j+1,k)-hy_starState(PRES_VAR,i,j-1,k)) &
!                   / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!             flatTilde(JAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!             if (hy_starState(VELY_VAR,i,j+1,k)<hy_starState(VELY_VAR,i,j,k)) then
!                flatTilde(JAXIS,i,j,k) = max(flatTilde(JAXIS,i,j,k), &
!                      &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!             else
!                flatTilde(JAXIS,i,j,k) = 1.0
!             end if
! #if NDIM==3
!             beta = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
!                   / max(1.e-99,abs(hy_starState(PRES_VAR,i,j,k+2)-hy_starState(PRES_VAR,i,j,k-2)))
!             Z    = abs(hy_starState(PRES_VAR,i,j,k+1)-hy_starState(PRES_VAR,i,j,k-1)) &
!                   / (hy_starState(GAMC_VAR,i,j,k)*hy_starState(PRES_VAR,i,j,k))
!             flatTilde(KAXIS,i,j,k) = max(0.,min(1.,(betaMax-beta)/(betaMax-betaMin)))
!             if (hy_starState(VELZ_VAR,i,j,k+1)<hy_starState(VELZ_VAR,i,j,k)) then
!                flatTilde(KAXIS,i,j,k) = max(flatTilde(KAXIS,i,j,k), &
!                      &                   min(1., (Zmax-Z)/(Zmax-Zmin)))
!             else
!                flatTilde(KAXIS,i,j,k) = 1.0
!             end if
! #endif
! #endif
!          end do
!       end do
!    end do
!    do k=limits(LOW,KAXIS)-kz, limits(HIGH,KAXIS)+kz
!       do j=limits(LOW,JAXIS)-ky, limits(HIGH,JAXIS)+ky
!          do i=limits(LOW,IAXIS)-kx, limits(HIGH,IAXIS)+kx
!             flat3d(i,j,k) = minval(flatTilde(1:NDIM,i,j,k))
! #ifdef FLAT_VAR
!             hy_starState(FLAT_VAR,i,j,k) = flat3d(i,j,k)
! #endif
!          end do
!       end do
!    end do
! end subroutine flattening
 
subroutine ensurePhysicalState(state)
   use Hydro_data, ONLY : hy_smalldens, hy_smallE, hy_smallpres, hy_smallX
   implicit none
   real, dimension(NRECON), intent(INOUT), target  :: state
   integer :: s
   real :: spcSumInv
   real, pointer :: spc(:)
   !$omp declare target
   state(HY_DENS) = max(hy_smalldens, state(HY_DENS))
   state(HY_PRES) = max(hy_smallpres, state(HY_PRES))
#if NSPECIES>0
   ! Limit and renormalize the species.
   spc => state(HY_NUM_VARS+1:HY_NUM_VARS+NSPECIES)
   do s = 1, NSPECIES
      spc(s) = max(hy_smallX,min(1.0,spc(s)))
   end do
   spcSumInv = 1./sum(spc(1:NSPECIES))
   spc = spc*spcSumInv
#endif
end subroutine ensurePhysicalState
 
 end subroutine hy_rk_getFaceFlux
 
 subroutine get_scratch_indices(dirLims, i_new, j_new, k_new, i, j, k, guardCells)
   implicit none
   integer ,intent(in) :: i, j ,k
   integer, dimension(3) :: guardCells
   integer ,intent(in), dimension(LOW:HIGH,MDIM)  ::  dirLims
   integer ,intent(OUT) :: i_new, j_new, k_new
   !$omp declare target
   i_new = 1 + i - dirLims(LOW,1) + guardCells(1)
   j_new = 1 + j - dirLims(LOW,2) + guardCells(2)
   k_new = 1 + k - dirLims(LOW,3) + guardCells(3)
 
 end subroutine get_scratch_indices

 subroutine initializeSnake(dir, dirLims, guardCells)
   use Hydro_data, only: snake,grv,shck,flat,hy_starState, uPlusArray, uMinusArray,flat3d,hy_grav
implicit none
integer :: dir, l,i,j,k,n
integer, dimension(3) :: guardCells
integer :: i_s, j_s, k_s   ! Scratch indices
integer, dimension(LOW:HIGH,MDIM), intent(IN) :: dirLims
!!$omp declare target

! Give snake proper values
! Maps the hk_solnData to a snake, such that it is contiguous in memory in the direction of iteration.
 !$omp target teams distribute parallel do collapse(3) default(none) shared(dir,hy_starState,snake,grv,hy_grav,shck,flat,flat3d,dirLims,guardCells) private(i,j,k,n,i_s, j_s, k_s)
   do k = dirLims(LOW,3) - guardCells(3), dirLims(HIGH,3) + guardCells(3)
      do j = dirLims(LOW,2) - guardCells(2), dirLims(HIGH,2) + guardCells(2)
         do i = dirLims(LOW,1) - guardCells(1), dirLims(HIGH,1) + guardCells(1)
            call get_scratch_indices(dirLims,i_s,j_s,k_s,i,j,k, guardCells)
            if (dir == IAXIS) then
            snake(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,i,j,k)
            snake(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,i,j,k)
            snake(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,i,j,k)
            snake(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,i,j,k)
            snake(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,i,j,k)
            snake(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,i,j,k)
            snake(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,i,j,k)*hy_starState(EINT_VAR,i,j,k)
#ifdef SPARK_GLM
            snake(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,i,j,k)
            snake(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,i,j,k)
            snake(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,i,j,k)
            snake(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,i,j,k)
#endif
#if NSPECIES+NMASS_SCALARS>0
            do n=SPECIES_BEGIN, MASS_SCALARS_END
               snake(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,i,j,k)
            enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
            grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,i,j,k)
#else
            grv(i_s,j_s,k_s) = hy_grav(IAXIS,i,j,k)
#endif
#endif
#ifdef SHOK_VAR
            shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,i,j,k)
#else
            shck(i_s,j_s,k_s) = 0.0
#endif
            flat(i_s,j_s,k_s) = flat3d(i,j,k)
            
else if (dir == JAXIS) then
            snake(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,i,k)
            snake(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,j,i,k)
            snake(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,j,i,k)
            snake(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,j,i,k)
            snake(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,j,i,k)
            snake(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,j,i,k)
            snake(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,i,k)*hy_starState(EINT_VAR,j,i,k)
#ifdef SPARK_GLM
            snake(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,j,i,k)
            snake(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,j,i,k)
            snake(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,j,i,k)
            snake(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,j,i,k)
#endif
#if NSPECIES+NMASS_SCALARS>0
            do n=SPECIES_BEGIN, MASS_SCALARS_END
               snake(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,j,i,k)
            enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
            grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,j,i,k)
#else
            grv(i_s,j_s,k_s) = hy_grav(IAXIS,j,i,k)
#endif
#endif
#ifdef SHOK_VAR
            shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,j,i,k)
#else
            shck(i_s,j_s,k_s) = 0.0
#endif
            flat(i_s,j_s,k_s) = flat3d(j,i,k)

else if (dir == KAXIS) then
            snake(HY_DENS,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,k,i)
            snake(HY_VELX,i_s,j_s,k_s) = hy_starState(VELX_VAR,j,k,i)
            snake(HY_VELY,i_s,j_s,k_s) = hy_starState(VELY_VAR,j,k,i)
            snake(HY_VELZ,i_s,j_s,k_s) = hy_starState(VELZ_VAR,j,k,i)
            snake(HY_PRES,i_s,j_s,k_s) = hy_starState(PRES_VAR,j,k,i)
            snake(HY_GAMC,i_s,j_s,k_s) = hy_starState(GAMC_VAR,j,k,i)
            snake(HY_RHOE,i_s,j_s,k_s) = hy_starState(DENS_VAR,j,k,i)*hy_starState(EINT_VAR,j,k,i)
#ifdef SPARK_GLM
            snake(HY_MAGX,i_s,j_s,k_s) = hy_starState(MAGX_VAR,j,k,i)
            snake(HY_MAGY,i_s,j_s,k_s) = hy_starState(MAGY_VAR,j,k,i)
            snake(HY_MAGZ,i_s,j_s,k_s) = hy_starState(MAGZ_VAR,j,k,i)
            snake(HY_PSIB,i_s,j_s,k_s) = hy_starState(PSIB_VAR,j,k,i)
#endif
#if NSPECIES+NMASS_SCALARS>0
            do n=SPECIES_BEGIN, MASS_SCALARS_END
               snake(HY_NUM_VARS+1+n-SPECIES_BEGIN,i_s,j_s,k_s)    = hy_starState(n,j,k,i)
            enddo
#endif
#ifdef GRAVITY
#ifdef GPOT_VAR
            grv(i_s,j_s,k_s) = hy_starState(GPOT_VAR,j,k,i)
#else
            grv(i_s,j_s,k_s) = hy_grav(IAXIS,j,k,i)
#endif
#endif
#ifdef SHOK_VAR
            shck(i_s,j_s,k_s) = hy_starState(SHOK_VAR,j,k,i)
#else
            shck(i_s,j_s,k_s) = 0.0
#endif
            flat(i_s,j_s,k_s) = flat3d(j,k,i)
         endif

         enddo
      enddo
   enddo



end subroutine initializeSnake


#include "reconstruct.F90"
#include "riemann.F90"
