!!****ih* source/physics/Hydro/HydroMain/unsplit/hy_uhdInterface
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
!!  hy_uhdInterface
!!
!! SYNOPSIS
!!  use hy_uhdInterface
!!
!! DESCRIPTION
!!  This Fortran module contains explicit interfaces for private
!!  subroutines specific to the "unsplit" solver(s).
!!  The latter correspond to the following Hydro "unsplit"
!!  implementations in FLASH4:
!!   1. MHD_StaggeredMesh
!!   2. Hydro_Unsplit
!!
!!  None of the subroutines below should be referred to directly
!!  by any code outside of the Hydro unit.
!!
!!***

Module hy_uhdInterface

  use hy_interface, ONLY : hy_setMinTimeStep, hy_shockDetect

  implicit none


#include "constants.h"
#include "Simulation.h"
#include "UHD.h"

  interface
     subroutine hy_avgState(sweepDir,VL,VR,Vavg)
       implicit none
       integer, intent(IN) :: sweepDir
       real, dimension(HY_VARINUM3), intent(IN)  :: VL,VR
       real, dimension(HY_VARINUM2), intent(OUT) :: Vavg
     end subroutine hy_avgState
  end interface




  interface
     subroutine hy_entropyFix(lambda,lambdaL,lambdaR)
       implicit none
       real, dimension(HY_WAVENUM), intent(INOUT) :: lambda
       real, dimension(HY_WAVENUM), intent(IN)    :: lambdaL,lambdaR
     end subroutine hy_entropyFix
  end interface



  interface
     Subroutine hy_dataReconstOneStep(tileDesc,U,loGC,hiGC,order,ix,iy,iz, &
                                          dt,del,ogravX,ogravY,ogravZ,&
                                          DivU,FlatCoeff,   &
                                          TransX_updateOnly,&
                                          TransY_updateOnly,&
                                          TransZ_updateOnly,&
                                          Wp,Wn,sig,&
                                          lambda,leftEig,rghtEig,&
                                          cellCfl,&
                                          hy_SpcR,hy_SpcL,hy_SpcSig)
       use Grid_tile,   ONLY : Grid_tile_t

       implicit none
       type(Grid_tile_t),intent(IN) :: tileDesc
       real,dimension(:,:,:,:),pointer :: U
       integer, intent(IN), dimension(MDIM):: loGC, hiGC
       integer,intent(IN) :: order,ix,iy,iz
       real,   intent(IN) :: dt
       real,   intent(IN), dimension(MDIM) :: del
       real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), &
            intent(IN), target :: ogravX,ogravY,ogravZ
       real, dimension(NDIM,loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), &
            intent(IN) :: FlatCoeff
       real, dimension(loGC(IAXIS):hiGC(IAXIS),loGC(JAXIS):hiGC(JAXIS),loGC(KAXIS):hiGC(KAXIS)), &
            intent(IN) :: DivU
       logical, intent(IN) ::  TransX_updateOnly, TransY_updateOnly, TransZ_updateOnly
       real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: Wp, Wn
       real, dimension(HY_VARINUMMAX,           NDIM),intent(OUT) :: sig
       real, dimension(              HY_WAVENUM,NDIM),intent(OUT) :: lambda
       real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: leftEig
       real, dimension(HY_VARINUM,   HY_WAVENUM,NDIM),intent(OUT) :: rghtEig
       real, optional, intent(INOUT) :: cellCfl ! We may pass back a lowered CFL factor here - KW
       real, pointer, optional, dimension(:,:,:,:,:) :: hy_SpcR,hy_SpcL,hy_SpcSig
     end subroutine hy_dataReconstOnestep
  end interface



  interface
     Subroutine hy_DataReconstructNormalDir_MH&
          (dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_DataReconstructNormalDir_MH
  end interface




  interface
     Subroutine hy_DataReconstructNormalDir_PPM&
          (dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_DataReconstructNormalDir_PPM
  end interface



  interface
     Subroutine hy_DataReconstructNormalDir_WENO&
          (wenoMethod,dir,dt,delta,Data1D,DataGrav1D,&
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: wenoMethod,dir
       real,   intent(IN) :: dt,delta
       real, pointer, dimension(:,:) :: Data1D
       real, pointer, dimension(:)   :: DataGrav1D
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     End Subroutine Hy_DataReconstructNormalDir_WENO
  end interface



  interface
     Subroutine hy_DataReconstructNormalDir_GP&
          (dir,dt,delta,DataMultiD,DataGravMultiD,&
          x1,x2,x3,radius, &
          FlatCoeff,TransUpdateOnly, &
          lambda0,leig0,reig0,&
          Wp,Wm,sig,&
          dnBnormal,aBnormal,& !These are optional
          dnGLMnormal,aGLMnormal,&
          Sr,Sl,SpcSig)        !These are optional
       implicit none
       integer,intent(IN) :: dir
       real,   intent(IN) :: dt,delta
#if NDIM < 3
       ! Although we define this array for 1D here for the sake of compilations,
       ! GP does not support 1D calculations.
       real, pointer, dimension(:,:,:)   :: DataMultiD,DataGravMultiD
#elif NDIM == 3
       real, pointer, dimension(:,:,:,:) :: DataMultiD,DataGravMultiD
#endif
       real, pointer, dimension(:) :: x1
       real, pointer, dimension(:) :: x2
       real, pointer, dimension(:) :: x3
       integer, intent(IN) :: radius
       real,    intent(IN) :: FlatCoeff
       logical, intent(IN) :: TransUpdateOnly
       real, dimension(           HY_WAVENUM),intent(OUT) :: lambda0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: leig0
       real, dimension(HY_VARINUM,HY_WAVENUM),intent(OUT) :: reig0
       real,intent(OUT),dimension(HY_VARINUMMAX) :: Wp,Wm
       real,intent(OUT),dimension(HY_VARINUMMAX) :: sig
       !optional arguments for MHD-----
       real,intent(IN), optional :: dnBnormal,dnGLMnormal
       real,intent(IN), dimension(HY_VARINUM), optional :: aBnormal,aGLMnormal
       !-------------------------------
       real,intent(OUT),dimension(HY_NSPEC),   optional :: Sr,Sl
       real,intent(OUT),dimension(:),          optional :: SpcSig
     end Subroutine hy_DataReconstructNormalDir_GP
    end interface



  interface
     subroutine hy_Roe(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm,Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_Roe
  end interface



  interface
     subroutine hy_HLL(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1),  intent(OUT) :: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_HLL
  end interface



  interface
     subroutine hy_HLLC(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_HLLC
  end interface



  interface
     subroutine hy_Marquina(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_Marquina
  end interface


  interface
     subroutine hy_MarquinaModified(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_MarquinaModified
  end interface


  interface
     subroutine hy_LLF(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_LLF
  end interface


!!$  interface
!!$     subroutine hy_upwindTransverseFlux(order,vector,lambda,leig,reig,sig)
!!$       implicit none
!!$       integer,intent(IN) :: order
!!$       real,intent(IN),dimension(HY_VARINUMMAX,-2:2)   :: vector
!!$       real,intent(IN),dimension(HY_WAVENUM) :: lambda
!!$       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
!!$       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
!!$       real,intent(OUT), dimension(HY_VARINUMMAX) :: sig
!!$     end subroutine hy_upwindTransverseFlux
!!$  end interface


  interface
     subroutine hy_upwindTransverseFlux&
          (dir,order,vm1,vc0,vp1,lambda,leig,reig,sigSize,sig,speciesScalar)
       implicit none
       integer,intent(IN) :: dir,order
       real,pointer,dimension(:)  :: vm1,vc0,vp1
       real,intent(IN),dimension(HY_WAVENUM) :: lambda
       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: leig
       real,intent(IN),dimension(HY_VARINUM,HY_WAVENUM) :: reig
!!$       real,pointer,dimension(:) :: lambda
!!$       real,pointer,dimension(:,:) :: leig,reig
       integer,intent(IN) :: sigSize
       real,intent(OUT),dimension(sigSize) :: sig
       logical, intent(IN),optional :: speciesScalar
     end subroutine hy_upwindTransverseFlux
  end interface



  interface
     subroutine hy_TVDslope(dir,VL,V0,VR,lambda,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX),intent(IN)  :: VL,V0,VR
       real, dimension(HY_WAVENUM),   intent(IN)  :: lambda
       real, dimension(HY_VARINUM,HY_WAVENUM), intent(IN) :: leig
       real, dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_TVDslope
  end interface


  interface
     subroutine hy_TVDslopeUpwind(dir,VLL,VL,V0,VR,VRR,lambdaL,lambda,lambdaR,leig,delbar)
       implicit none
       integer, intent(IN) :: dir
       real,dimension(HY_VARINUMMAX),intent(IN)  :: VLL,VL,V0,VR,VRR
       real,dimension(HY_WAVENUM),  intent(IN)  :: lambdaL,lambda,lambdaR
       real,dimension(HY_VARINUM,HY_WAVENUM),intent(IN) :: leig
       real,dimension(HY_VARINUMMAX),intent(OUT) :: delbar
     end subroutine hy_TVDslopeUpwind
  end interface




  interface
     subroutine hy_addViscousFluxes&
          (tileDesc,blkLimitsGC,ix,iy,iz,Flux,mu,sweepDir)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)   :: tileDesc
       integer, INTENT(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: mu
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                  intent(IN) :: mu
#endif
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_addViscousFluxes
    end interface



  interface
     subroutine hy_addThermalFluxes&
          (tileDesc,blkLimitsGC,ix,iy,iz,Flux,kappa,sweepDir)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)   :: tileDesc
       integer, INTENT(IN) :: ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: kappa
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                  intent(IN) :: kappa
#endif
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_addThermalFluxes
    end interface



    interface
       subroutine hy_eigenParameters&
            (V,dir,U_normal,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)
         implicit none
         real, dimension(HY_VARINUM2), intent(IN)  :: V
         integer, intent(IN)  :: dir
         real, intent(OUT) :: U_normal,C_fast
         real, intent(OUT), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(OUT),optional :: B_beta
         real, intent(OUT), optional :: C_hyp
       end subroutine hy_eigenParameters
    end interface



    interface
       subroutine hy_eigenValue(EigValue,U_normal,C_fast,C_alfn,C_slow,C_hyp)
         implicit none
         real,dimension(HY_WAVENUM), intent(OUT) :: EigValue
         real,intent(IN) :: U_normal,C_fast
         real,intent(IN), optional :: C_alfn,C_slow,C_hyp
       end subroutine hy_eigenValue
    end interface



    interface
       subroutine hy_eigenVector&
            (LeftEigvec,RightEigvec,V,dir,cons,C_fast,C_alfn,C_slow,A_f,A_s,B_beta,C_hyp)
         implicit none
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: LeftEigvec
         real, dimension(HY_VARINUM,HY_WAVENUM), intent(OUT) :: RightEigvec
         real, dimension(HY_VARINUM2), intent(IN) :: V
         integer, intent(IN) :: dir
         logical, intent(IN) :: cons
         real,    intent(IN) :: C_fast
         real,    intent(IN), optional :: C_alfn,C_slow,A_f,A_s
         real, dimension(MDIM), intent(IN), optional  :: B_beta
         real,    intent(IN), optional :: C_hyp
       end subroutine hy_eigenVector
    end interface



    interface
       subroutine hy_prim2con(V,CU)
         implicit none
         real ,dimension(HY_VARINUM2), intent(IN) :: V
         real ,dimension(HY_VARINUM),  intent(OUT) :: CU
       end subroutine hy_prim2con
    end interface



    interface
       subroutine hy_con2prim(CU,game,V)
         implicit none
         real ,dimension(HY_VARINUM), intent(IN)  :: CU
         real, intent(IN) :: game
         real ,dimension(HY_VARINUM), intent(OUT) :: V
       end subroutine hy_con2prim
    end interface



    interface
       subroutine hy_prim2flx(dir,V,F)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(*),            intent(IN) :: V
         real, dimension(HY_VARINUM1),  intent(OUT) :: F
       end subroutine hy_prim2flx
    end interface



    interface
       Subroutine hy_checkRHjumpCond(dir,dens,velocity,pres,gamc,Wp,Wn,SWp,SWn)
         implicit none
         integer, intent(IN)  :: dir
         real, dimension(MDIM), intent(IN) :: velocity
         real, intent(IN) :: dens, pres, gamc
         real, dimension(HY_DENS:HY_EINT), intent(IN) :: Wp,Wn
         logical, intent(OUT) :: SWp,SWn
       end Subroutine hy_checkRHjumpCond
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


!! FOR UNSPLIT STAGGERED MESH MHD SOLVER -------------------------------------------
#if defined(FLASH_USM_MHD) || defined(FLASH_UGLM_MHD)
  interface
     subroutine hy_addResistiveFluxes&
          (tileDesc,blkLimitsGC,ix,iy,iz,Flux,eta,sweepDir)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       type(Grid_tile_t), intent(IN)   :: tileDesc
       integer, INTENT(IN) :: block,ix,iy,iz
       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
#ifdef FIXEDBLOCKSIZE
  real, dimension(GRID_ILO_GC:GRID_IHI_GC, &
                  GRID_JLO_GC:GRID_JHI_GC, &
                  GRID_KLO_GC:GRID_KHI_GC),intent(IN) :: eta
#else
  real, dimension(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS), &
                  blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS), &
                  blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)),&
                  intent(IN) :: eta
#endif
       integer, INTENT(IN) :: sweepDir
     end subroutine hy_addResistiveFluxes
    end interface





  interface
     subroutine hy_HLLD(dir,Vm,Vp,Fstar,speed,ierr)
       implicit none
       integer, intent(IN) :: dir
       real, dimension(HY_VARINUMMAX), intent(IN) :: Vm, Vp
       real, dimension(HY_VARINUM1), intent(OUT):: Fstar
       real, intent(OUT) :: speed
       integer, intent(OUT) :: ierr
     end subroutine hy_HLLD
  end interface





  interface
     Subroutine hy_addBiermannBatteryTerms(tileDesc,blkLimitsGC,ix,iy,iz,Flux,sweepDir)
       use Grid_tile, ONLY : Grid_tile_t
       implicit none
       integer, INTENT(IN) :: ix,iy,iz
       type(Grid_tile_t), intent(IN) :: tileDesc

       integer, dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimitsGC
       real, dimension(HY_VARINUM), intent(INOUT) :: Flux
       integer, INTENT(IN) :: sweepDir
     end Subroutine hy_addBiermannBatteryTerms
  end interface

#endif
!endif #ifdef FLASH_USM_MHD


End Module hy_uhdInterface
