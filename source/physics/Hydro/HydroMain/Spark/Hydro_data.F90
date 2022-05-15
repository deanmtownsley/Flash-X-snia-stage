!!****if* source/physics/Hydro/HydroMain/Spark/Hydro_data
!!
!!  NAME
!!    Hydro_data
!!
!!  SYNOPSIS
!!    use Hydro_data
!!
!!  DESCRIPTION
!!    Stores data for Spark Hydro
!!
!!  NOTES
!!
!!***
!!REORDER(4): hy_fluxBuf[XYZ]

module Hydro_data
#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  implicit none
  save
  !$omp declare target (dirLims,flat3d,hy_starState,uPlusArray,uMinusArray,snake,solState_tmp)
  !$omp declare target (flat,grv,shck,xsize,ysize,zsize,hy_tiny,hy_hybridRiemann,hy_geometry)
  !$omp declare target (hy_C_hyp,flux,hy_smalldens, hy_smallE, hy_smallpres, hy_smallX,hy_cvisc,del)
  !$omp declare target (hy_flx,hy_fly,hy_flz,hy_alphaGLM,hy_grav)

  ! Added for GPU
  real, allocatable, target :: flux(:,:,:,:)
  real, save, allocatable, target :: snake(:,:,:,:)
  real, save, target, allocatable :: uPlusArray(:,:,:,:)
  real, save, target, allocatable :: uMinusArray(:,:,:,:)
  real, save, allocatable :: flat(:,:,:)
  real, save, allocatable :: flat3d(:,:,:)
  integer :: xsize,ysize,zsize
  real, save, allocatable :: grv(:,:,:)
  real, save, allocatable :: shck(:,:,:)
  integer, save, dimension(LOW:HIGH,MDIM) :: dirLims
  logical :: hydro_GPU_scratch = .False.
  real, dimension(MDIM) :: del
  real, allocatable, target :: solState_tmp(:,:,:,:)
  logical :: scratch_allocated


  integer :: hy_meshMe
  real :: hy_cfl
  logical :: hy_hydroComputeDtFirstCall
  logical :: hy_updateHydroFluxes
  real :: hy_dt, hy_dtmin
  integer :: hy_gcMaskSize
  integer :: hy_globalComm
  integer :: hy_meshNumProcs
  logical :: hy_restart
  logical :: hy_shockDetectOn
  logical :: hy_useTiling
  real :: hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_smallu

  !One block's worth of fluxes defined generally (not assuming fixed block size mode)
  real, allocatable, target :: hy_flx(:,:,:,:), hy_fly(:,:,:,:), hy_flz(:,:,:,:)
  
  !Flux buffers
#ifdef FIXEDBLOCKSIZE
  real, dimension(NFLUXES, GRID_ILO:GRID_IHI+1  , &
    & GRID_JLO:GRID_JHI    ,GRID_KLO:GRID_KHI  ),target :: hy_fluxBufX
# if NDIM > 1
  real, dimension(NFLUXES, GRID_ILO:GRID_IHI    , &
    & GRID_JLO:GRID_JHI+1  ,GRID_KLO:GRID_KHI  ),target :: hy_fluxBufY
# else
  real, dimension(NFLUXES, 0    ,0  ,0  ) :: hy_fluxBufY
# endif
# if NDIM > 2
  real, dimension(NFLUXES, GRID_ILO:GRID_IHI    , &
    & GRID_JLO:GRID_JHI    ,GRID_KLO:GRID_KHI+1),target :: hy_fluxBufZ
# else
  real, dimension(NFLUXES, 0    ,0  ,0  ) :: hy_fluxBufZ
# endif
# else  
  !Not fixed block size
  real, allocatable, target :: hy_fluxBufX(:,:,:,:),hy_fluxBufY(:,:,:,:),hy_fluxBufZ(:,:,:,:)
# endif
 
  real, allocatable :: hy_grav(:,:,:,:)  

  logical :: hy_useHydro
  logical :: hy_fluxCorrect, hy_fluxCorrectPerLevel
  integer, dimension(NFLUXES) :: hy_fluxCorVars
  integer :: hy_geometry

  logical :: hy_threadWithinBlock
  logical, dimension(NUNK_VARS) :: hy_gcMask
  ! Additional scratch storage for RK time stepping
  real, allocatable, target :: hy_starState(:,:,:,:)

  ! Limiter info
  real :: hy_limRad
  real :: hy_cvisc

  real :: hy_tiny=1.e-32
  real :: hy_gravConst, hy_4piGinv

  logical :: hy_hybridRiemann, hy_flattening

  real :: hy_C_hyp, hy_alphaGLM, hy_lChyp
  real :: hy_bref
  ! System of units used
  character(4) :: hy_units

end module Hydro_data
