!!****ih* source/physics/Hydro/localAPI/hy_interface
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
!!  hy_interface
!!
!! SYNOPSIS
!!  use hy_interface
!!
!! DESCRIPTION
!!  This Fortran module contains explicit interfaces for private
!!  subroutines in various implementations of the HydroMain subunit.
!!
!!  These interfaces are shared within the Hydro unit: subroutines
!!  implemented in one subtree of the Hydro unit could be referred
!!  to by code in a different subtree.
!!
!!  None of the subroutines below should be referred to directly
!!  by any code outside of the Hydro unit.
!!
!!***

!!REORDER(4): [xyz]flux, faceFlux, fl[xyz], fluxBuf[XYZ], Sp, U

Module hy_interface

  implicit none


#include "constants.h"
#include "Simulation.h"

  interface
     subroutine hy_hllComputeFluxes ( tileLimits, Uin, plo, flX, flY, flZ, loFl, del, dt )
       implicit none
       integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
       integer, intent(IN)  :: plo(*)
       real,    intent(IN)  :: Uin(plo(1):,plo(2):,plo(3):,plo(4):)
       integer, intent(IN)  :: loFl(*)
       real,    intent(OUT) :: flX(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
       real,    intent(OUT) :: flY(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
       real,    intent(OUT) :: flZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
       real,    intent(IN)  :: del(1:MDIM)
       real,    intent(IN)  :: dt
     end subroutine hy_hllComputeFluxes
  end interface

  interface
     subroutine hy_hllUpdateSolution( tileLimits, Uin, plo, Uout, flX, flY, flZ, loFl, del, dt )
       implicit none
       integer, intent(IN)  :: tileLimits(LOW:HIGH, 1:MDIM)
       integer, intent(IN)  :: plo(*)
       real,    intent(INOUT):: Uin(plo(1):,plo(2):,plo(3):,plo(4):)
       real,    intent(OUT) :: Uout(plo(1):,plo(2):,plo(3):,plo(4):)
       integer, intent(IN)  :: loFl(*)
       real,    intent(IN)  :: flX(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
       real,    intent(IN)  :: flY(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
       real,    intent(IN)  :: flZ(loFl(1):,loFl(2):,loFl(3):,loFl(4):)
       real,    intent(IN)  :: del(1:MDIM)
       real,    intent(IN)  :: dt
       TARGET :: Uin,Uout
     end subroutine hy_hllUpdateSolution
  end interface



  interface
     subroutine hy_getRiemannState(tileDesc,U,blkLimits,loGC,hiGC,dt,del, &
                                       ogravX,ogravY,ogravZ,&
                                       scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,&
                                       hy_SpcR,&
                                       hy_SpcL,hy_SpcSig,normalFieldUpdate)

       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)   :: tileDesc
       integer, intent(IN),dimension(LOW:HIGH,MDIM):: blkLimits !, blkLimitsGC
       integer, intent(IN), dimension(MDIM):: loGC, hiGC
       real,    intent(IN)   :: dt
       real,    intent(IN),dimension(MDIM) :: del
!!$       real, dimension(blkLimitsGC(HIGH,IAXIS),  &
!!$                       blkLimitsGC(HIGH,JAXIS),  &
!!$                       blkLimitsGC(HIGH,KAXIS)), &
!!$                       intent(IN) :: ogravX,ogravY,ogravZ
  real, dimension(loGC(IAXIS):,loGC(JAXIS):,loGC(KAXIS):), intent(IN) :: ogravX,ogravY,ogravZ
       real,pointer,dimension(:,:,:,:)::U
       real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
       logical, intent(IN), optional :: normalFieldUpdate
     end subroutine hy_getRiemannState
 end interface



  interface
     subroutine hy_getFaceFlux (tileDesc,blkLimits,blkLimitsGC,del,&
                                     loFl, xflux,yflux,zflux,&
                                     scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr,scrch_Ptr,&
                                     hy_SpcR,hy_SpcL,hy_SpcSig,lastCall)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)  :: tileDesc
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
       real,    dimension(MDIM), intent(IN)         :: del
       integer, dimension(MDIM+1),intent(IN)        :: loFl
       real, intent(OUT) :: XFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZED TO PREVENT INDEX REORDERING!
       real, intent(OUT) :: YFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
       real, intent(OUT) :: ZFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
       real, pointer, dimension(:,:,:,:) :: scrchFaceXPtr,scrchFaceYPtr,scrchFaceZPtr
       real, pointer, dimension(:,:,:,:) :: scrch_Ptr
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
       logical, optional, intent(IN) :: lastCall
     end subroutine hy_getFaceFlux
  end interface



  interface
     subroutine hy_setMinTimeStep(tileDesc,i,j,k,delta,speed)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)   :: tileDesc
       integer, INTENT(in) :: i,j,k
       real, INTENT(in) :: delta,speed
     end subroutine hy_setMinTimeStep
  end interface




  interface
     subroutine hy_unsplitUpdate(tileDesc,Uin,Uout,rangeSwitch,dt,del,blkLimits,&    
                                     blGC,loFl,xflux,yflux,zflux,gravX,gravY,gravZ,&
                                     scrch_Ptr)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)   :: tileDesc
       real,dimension(:,:,:,:),pointer :: Uin, Uout
       integer,intent(IN) :: rangeSwitch
       real, intent(IN)   :: dt
       real, intent(IN)   :: del(MDIM)
       integer,intent(IN) :: blkLimits(LOW:HIGH,MDIM)
       integer,intent(IN) :: blGC(LOW:HIGH,MDIM)
       integer, dimension(MDIM+1),intent(IN)        :: loFl
       real, intent(IN) :: XFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZED TO PREVENT INDEX REORDERING!
       real, intent(IN) :: YFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
       real, intent(IN) :: ZFLUX (loFl(1):, loFl(2): ,loFl(3):, loFl(4): )!CAPITALIZATION INTENTIONAL!
    real, dimension(blGC(LOW,IAXIS):blGC(HIGH,IAXIS),blGC(LOW,JAXIS):blGC(HIGH,JAXIS),blGC(LOW,KAXIS):blGC(HIGH,KAXIS)),& 
         intent(IN) :: gravX,gravY,gravZ
       real, pointer, dimension(:,:,:,:) :: scrch_Ptr
     end subroutine hy_unsplitUpdate
  end interface

    interface
       subroutine hy_multiTempAfter(blockCount, blockList, dt)
         implicit none
         integer, intent(in) :: blockCount
         integer, intent(in) :: blockList(blockCount)
         real,    intent(in) :: dt
       end subroutine hy_multiTempAfter
    end interface

    interface
       subroutine hy_ragelike(uextra, soln, dens_old, &
            duele_adv, duion_adv, durad_adv, &
            xc, yc, zc, &
            uele_new, uion_new, urad_new, &
            stepFactor)
         implicit none
         real, intent(in)  :: uextra
         real, intent(in)  :: soln(NUNK_VARS)
         real, intent(in)  :: dens_old
         real, intent(in)  :: duele_adv
         real, intent(in)  :: duion_adv
         real, intent(in)  :: durad_adv
         real, intent(in)  :: xc 
         real, intent(in)  :: yc 
         real, intent(in)  :: zc
         real, intent(out) :: uele_new
         real, intent(out) :: uion_new
         real, intent(out) :: urad_new
         real, intent(in),OPTIONAL  :: stepFactor
       end subroutine hy_ragelike
    end interface

  interface
     subroutine hy_updateSpeciesMassScalar&
          (order,densNew,Sp,U,FL,FR,GL,GR,HL,HR,dx,dy,dz,dt,SpNew)
  implicit none
  integer, intent(IN) :: order
  real, intent(IN) :: densNew
#if NDIM == 1
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,1,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,1,1) :: U
#elif NDIM == 2
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,1)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,1) :: U
#elif NDIM == 3
  real, intent(IN), dimension(NSPECIES+NMASS_SCALARS,-3:3,-3:3,-3:3)   :: Sp
  real, intent(IN), dimension(6,-3:3,-3:3,-3:3)   :: U
#endif
  real, intent(IN)  :: FL,FR,GL,GR,HL,HR,dx,dy,dz,dt
  real, intent(OUT), dimension(NSPECIES+NMASS_SCALARS) :: SpNew
     end subroutine hy_updateSpeciesMassScalar
  end interface


  interface
     subroutine hy_energyFix(tileDesc,U,blkLimits,dt,dtOld,del,eosMode)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       
       type(Grid_tile_t), intent(IN)   :: tileDesc
       real, pointer, dimension(:,:,:,:) :: U
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
       real, intent(IN) :: dt,dtOld
       real, dimension(MDIM), intent(IN) :: del
       integer, intent(IN) :: eosMode
     end subroutine hy_energyFix
  end interface



  interface 
     subroutine hy_unitConvert(Uin,blkLimitsGC,convertDir)
       implicit none
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       integer, intent(IN) :: convertDir
       real, dimension(:,:,:,:) :: Uin
     end subroutine hy_unitConvert
  end interface
  


    interface
       subroutine hy_prepareNewGravityAccel(gcMaskLogged)
         logical,intent(in) :: gcMaskLogged
       end subroutine hy_prepareNewGravityAccel
    end interface


    interface
       subroutine hy_putGravity&
            (tileDesc,blGC,Uin,dt,dtOld,gravX,gravY,gravZ,potentialIndex,&
             lastCall)
         use Grid_tile, ONLY : Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(IN)   :: tileDesc
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blGC
         real,dimension(:,:,:,:),pointer :: Uin
         real,    intent(IN) :: dt, dtOld
         real,dimension(blGC(LOW,IAXIS):blGC(HIGH,IAXIS), blGC(LOW,JAXIS):blGC(HIGH,JAXIS), blGC(LOW,KAXIS):blGC(HIGH,KAXIS)),&
              intent(OUT) :: gravX,gravY,gravZ
         integer, intent(IN), OPTIONAL :: potentialIndex
         logical, intent(IN), OPTIONAL :: lastCall
       end subroutine hy_putGravity
    end interface



    interface
       Subroutine hy_addGravity&
            (tileDesc,blkLimits,loGC,hiGC,dt,gravX,gravY,gravZ)
         use Grid_tile, ONLY : Grid_tile_t
         implicit none
         type(Grid_tile_t), intent(IN)   :: tileDesc
         integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits
         integer, intent(IN), dimension(MDIM):: loGC, hiGC
         real,    intent(IN) :: dt
         real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), intent(IN) :: &
              gravX,gravY,gravZ
       end Subroutine hy_addGravity
    end interface



    interface
       subroutine hy_shockDetectBlk(Uin,loI,blkLimitsGC,Uout,loO,blkLimits,del )
       implicit none
       integer,dimension(LOW:HIGH,MDIM),INTENT(IN) :: blkLimits,blkLimitsGC
       integer, intent(IN)  :: loI(*), loO(*)
       real,    intent(IN)  :: UIN(loI(1):,loI(2):,loI(3):,loI(4):)
       real,    intent(OUT) :: UOUT(loO(1):,loO(2):,loO(3):,loO(4):)
!       real, dimension(:,:,:,:),INTENT(IN) :: Uin
!       real,dimension(:,:,:,:),INTENT(INOUT) :: Uout
       real,dimension(MDIM),INTENT(IN) :: del
     end subroutine hy_shockDetectBlk
  end interface

  interface
     subroutine hy_shockDetect()
       implicit none
     end subroutine hy_shockDetect
  end interface

    interface
     subroutine hy_gravityStep(simTime, dt, dtOld)
       implicit none
       real, intent(IN) ::  simTime, dt, dtOld
     end subroutine Hy_gravityStep
     
     subroutine hy_gravityStepBlk(tileDesc, blkLimitsGC, Uin, blkLimits, Uout, del,timeEndAdv,dt,dtOld)
       use Grid_tile, ONLY : Grid_tile_t
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       
       real, pointer, dimension(:,:,:,:) :: Uout
       real, pointer, dimension(:,:,:,:) :: Uin
       
       real,dimension(MDIM),intent(IN) :: del
       integer,dimension(LOW:HIGH,MDIM),intent(INoUt) ::blkLimits,blkLimitsGC 
       type(Grid_tile_t), intent(IN) :: tileDesc
     end subroutine hy_gravityStepBlk
  end interface

  ! Dongwook thinks that the following two interfaces should be taken care of for
  ! general implementation later.
!!$  interface
!!$     Subroutine hy_counterGP_init(radiusGP, counterGP, blkLimitsGC)
!!$       implicit none
!!$       real,    intent(IN)  :: radiusGP
!!$       integer, intent(OUT) :: counterGP
!!$       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
!!$     end subroutine hy_counterGP_init
!!$  end interface
!!$
  interface
     Subroutine hy_initGP(RinvGP, WpGP, WmGP, blkLimitsGC)
       implicit none
       real, intent(INOUT), dimension(:,:) :: RinvGP
       real, intent(INOUT), dimension(:,:) :: WpGP
       real, intent(INOUT), dimension(:,:) :: WmGP
       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimitsGC
     end Subroutine hy_initGP
  end interface

  interface hy_computeFluxes
     subroutine hy_computeFluxes(tileDesc, Uin, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder )
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t),intent(in) :: tileDesc
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
     end subroutine hy_computeFluxes
     subroutine hy_computeFluxes_fluxbuf(tileDesc, fluxBufX,fluxBufY,fluxBufZ,lo, Uin, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder )
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t),intent(in) :: tileDesc
       integer,intent(in) :: lo(3)
       real,intent(OUT),dimension(:, lo(1): ,lo(2): ,lo(3): ) :: fluxBufX,fluxBufY,fluxBufZ
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
     end subroutine hy_computeFluxes_fluxbuf
  end interface

    interface hy_updateSolution
     subroutine hy_updateSolution(tileDesc, Uin, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder )
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t),intent(in) :: tileDesc
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
     end subroutine hy_updateSolution

     subroutine hy_updateSolution_fluxbuf(tileDesc, fluxBufX,fluxBufY,fluxBufZ,lo, Uin, Uout, del,timeEndAdv, dt, dtOld,  &
          sweepOrder, updateMode)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t),intent(in) :: tileDesc
       integer,intent(in) :: lo(3)
       real,intent(IN),dimension(:, lo(1): ,lo(2): ,lo(3): ) :: fluxBufX,fluxBufY,fluxBufZ
       real, pointer, dimension(:,:,:,:) :: Uout,Uin
       real,    INTENT(IN) :: timeEndAdv, dt, dtOld
       integer, INTENT(IN) :: sweepOrder
       real,dimension(MDIM),intent(IN) :: del
       integer,intent(IN),OPTIONAL :: updateMode ! could be set to one of UPDATE_ALL, UPDATE_INTERIOR, UPDATE_BOUND
     end subroutine hy_updateSolution_fluxbuf
  end interface


!! FOR UNSPLIT STAGGERED MESH MHD SOLVER -------------------------------------------


  interface
     subroutine hy_addOhmicHeating(tileDesc,blkLimits,ix,iy,iz,Qohm,eta)
       use Grid_tile, ONLY : Grid_tile_t
     implicit none
     type(Grid_tile_t), intent(IN) :: tileDesc
     integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits
     integer, INTENT(IN) :: ix,iy,iz
     real, intent(INOUT) :: Qohm
     real, intent(IN)    :: eta  
     end subroutine hy_addOhmicHeating
  end interface


  interface
     subroutine hy_staggeredDivb(tileDesc,dt,del,blkLimits,blkLimitsGC,halfTimeAdvance)
       use Grid_tile, ONLY : Grid_tile_t 
       implicit none
       type(Grid_tile_t), intent(IN) :: tileDesc
       real,    intent(IN) :: dt
       real,    dimension(MDIM),   intent(IN) :: del
       integer, dimension(LOW:HIGH,MDIM), intent(IN) :: blkLimits,blkLimitsGC
       logical, intent(IN) :: halfTimeAdvance
     end subroutine hy_staggeredDivb
  end interface





  interface
     subroutine hy_getElectricFields(tileDesc,blkLimits,blkLimitsGC,del,flx,fly,flz)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN) :: tileDesc
       integer, intent(IN), dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
       real,    intent(IN), dimension(MDIM)  :: del  
#ifdef FIXEDBLOCKSIZE
       real, DIMENSION(NFLUXES,               &
                       GRID_ILO_GC:GRID_IHI_GC,  &
                       GRID_JLO_GC:GRID_JHI_GC,  &
                       GRID_KLO_GC:GRID_KHI_GC), &
                       intent(IN) :: flx,fly,flz
#else
       real, DIMENSION(NFLUXES,               &
                       blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                       blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                       blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                       intent(IN) :: flx,fly,flz
#endif
     end subroutine hy_getElectricFields
  end interface



  interface
     subroutine hy_getCurrents(tileDesc, rangeSwitch, blkLimits,datasize, del, Jp, Jm, &
                               mode_switch, scrch_Ptr,ix,iy,iz)
       use Grid_tile, ONLY : Grid_tile_t
        implicit none
         type(Grid_tile_t), intent(IN) :: tileDesc
        integer,intent(in) :: rangeSwitch
        integer,intent(in) :: blkLimits  (LOW:HIGH,MDIM)
        integer,intent(in) :: datasize(MDIM), mode_switch
        real,   intent(in) :: del(MDIM)
#ifdef FIXEDBLOCKSIZE
        real, intent(inout) :: Jp(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
        real, intent(inout) :: Jm(3,GRID_IHI_GC,GRID_JHI_GC,GRID_KHI_GC)
#else
        real, intent(inout) :: Jp(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))  
        real, intent(inout) :: Jm(3,datasize(IAXIS),datasize(JAXIS),datasize(KAXIS))
#endif
        real, pointer, dimension(:,:,:,:) :: scrch_Ptr
        integer, intent(in), OPTIONAL :: ix,iy,iz 
     end subroutine hy_getCurrents
  end interface



  interface
     subroutine hy_getFluxDeriv( ix,iy,iz,blkLimitsGC,&
                                     fluxType,DerivDir,   &
                                     faceFlux,            &
                                     Flux1Deriv,          &
                                     Flux2Deriv           )
       implicit none
       integer, intent(IN) :: ix,iy,iz
       integer, intent(IN) :: fluxType,DerivDir
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
#ifdef FIXEDBLOCKSIZE
       real,    dimension(NFLUXES,                  &
                          GRID_ILO_GC:GRID_IHI_GC,  &
                          GRID_JLO_GC:GRID_JHI_GC,  &
                          GRID_KLO_GC:GRID_KHI_GC), &
                          intent(IN) :: faceFlux
#else
       real,    dimension(NFLUXES,                  &
                          blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),  &
                          blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),  &
                          blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)), &
                          intent(IN) :: faceFlux
#endif
       real,    intent(OUT):: Flux1Deriv,Flux2Deriv
     end subroutine hy_getFluxDeriv
  end interface



  interface
     subroutine hy_biermannSource(blockCount, blockList, dt)
       implicit none
       integer, INTENT(IN) :: blockCount  
       integer, INTENT(IN), dimension(blockCount) :: blockList
       real,    INTENT(IN) :: dt
     end subroutine hy_biermannSource
  end interface

End Module hy_interface
