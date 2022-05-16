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
!!REORDER(4): hy_afluxBuf[XYZ]

module Hydro_data
#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  implicit none
  save
  !$omp declare target (hy_sizex,hy_sizey,hy_sizez,hy_tiny,hy_hybridRiemann,hy_geometry,hy_alphaGLM)
  !$omp declare target (hy_C_hyp,hy_smalldens, hy_smallE, hy_smallpres, hy_smallX,hy_cvisc,hy_del)
  !$omp declare target (hy_ablk,hy_autmphy_aflx,hy_afly,hy_aflz,hy_aflxbx,hy_aflxby,hy_aflxbz,hy_agrav, hy_aflat3d)
  !$omp declare target (hy_flux,hy_uplus, hy_uminus, hy_grv, hy_shk, hy_flat,hy_tposedBlk)
  !$omp declare target (hy_axcenter, hy_axright, hy_axleft, hy_aycenter, hy_azcenter, hy_aarea, hy_avol) 

  real, allocatable, dimension(:), target :: hy_ablk,hy_atposedBlk,hy_aflx,hy_afly,&
       hy_aflz,hy_aflxbx,hy_aflxby,hy_aflxbz, hy_agrav, hy_aflat3d, hy_autmp
  real,  allocatable, dimension(:), target :: hy_axcenter, hy_axright, hy_axleft, hy_aycenter, hy_azcenter, hy_aarea, hy_avol 
  real, allocatable, dimension(:,:,:,:), target :: hy_flux, hy_flat,  hy_uplus, hy_uminus, hy_grv, hy_shk,hy_transBlk
  real,  pointer, dimension(:,:,:,:) :: hy_starState, hy_tmpState
!!$  real, pointer, dimension(:,:,:,:) hy_uplus, hy_uminus
  real,  pointer, dimension(:,:,:,:) :: hy_flx,hy_fly,hy_flz,hy_flxbx,hy_flxby,hy_flxbz,hy_flux
  real, pointer, dimension(:,:,:) :: hy_area, hy_vol, hy_flat3d
!!$  real,  pointer, dimension(:,:,:) :: hy_flat, hy_grv, hy_grav,hy_shk, hy_tposedBlk
  real,  pointer, dimension(:) :: hy_xcenter, hy_ycenter, hy_zcenter, hy_xright, hy_xleft


  !! These are the variable that are offloaded
  integer :: hy_sizex,hy_sizey,hy_sizez
  integer :: hy_geometry
  integer, dimension(LOW:HIGH,MDIM) :: dirLims

  real :: hy_tiny=1.e-32
  real :: hy_C_hyp, hy_alphaGLM, hy_lChyp
  real :: hy_bref
  real :: hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_smallu
  real :: hy_cvisc
  real :: hy_del(MDIM)

  logical :: hy_hybridRiemann, hy_flattening
  

  !! AD: These variables don't see to be
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
  logical :: hy_useHydro
  logical :: hy_fluxCorrect, hy_fluxCorrectPerLevel
  integer, dimension(NFLUXES) :: hy_fluxCorVars
  real :: hy_limRad

  real :: hy_gravConst, hy_4piGinv



  logical :: hy_threadWithinBlock
  logical, dimension(NUNK_VARS) :: hy_gcMask

  
end module Hydro_data
