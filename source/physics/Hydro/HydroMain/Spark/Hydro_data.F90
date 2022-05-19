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


module Hydro_data
#include "constants.h"
#include "Simulation.h"
#include "Spark.h"

  implicit none
  save
  !$omp declare target (hy_sizex,hy_sizey,hy_sizez,hy_tiny,hy_hybridRiemann,hy_geometry,hy_alphaGLM)
  !$omp declare target (hy_C_hyp,hy_smalldens, hy_smallE, hy_smallpres, hy_smallX,hy_cvisc,hy_del)
  !$omp declare target (hya_xcenter, hya_xright, hya_xleft, hya_ycenter, hya_zcenter, hya_area, hya_vol) 
  !$omp declare target (hy_xcenter, hy_xright, hy_xleft, hy_ycenter, hy_zcenter, hy_area, hy_vol) 
  !$omp declare target (hya_starState, hya_tmpState,hya_flx,hya_fly)
  !$omp declare target (hya_flz,hya_flxbx,hya_flxby,hya_flxbz, hya_grav, hya_utmp, hya_flat3d)
  !$omp declare target (hy_starState,hy_tmpState,hy_flx,hy_fly,hy_flz,hy_flxbx,hy_flxby,hy_flxbz, hy_grav,hy_flat3d)
  !$omp declare target (hy_flux,  hy_uplus, hy_uminus,hy_tposeBlk, hy_grv, hy_shk, hy_flat)

  !! Arrays and their pointers that are needed only for non cartesian geometry
  real,  allocatable, dimension(:), target :: hya_xcenter, hya_xright, hya_xleft, &
       hya_ycenter, hya_zcenter, hya_area, hya_vol 
  real,  pointer, dimension(:) :: hy_xcenter, hy_ycenter, hy_zcenter, hy_xright, hy_xleft
  real, pointer, dimension(:,:,:) :: hy_area, hy_vol

  !! Arrays and their pointers that will always have global indices
  real, allocatable, dimension(:), target :: hya_starState, hya_tmpState,hya_flx,hya_fly,&
       hya_flz,hya_flxbx,hya_flxby,hya_flxbz, hya_grav, hya_flat3d
  real, pointer, dimension(:,:,:,:) :: hy_starState,hy_tmpState,hy_flx,hy_fly,&
       hy_flz,hy_flxbx,hy_flxby,hy_flxbz, hy_grav
  real, pointer,dimension(:,:,:) :: hy_flat3d

  !! Arrays that are beginning life with local indices for now, with their pointers commented out
  real, allocatable, dimension(:,:,:,:), target :: hy_flux,  hy_uplus, hy_uminus,hy_tposeBlk
  real, allocatable, dimension(:,:,:), target :: hy_grv, hy_shk, hy_flat
!!$  real, allocatable, dimension(:), target :: hya_flux,  hya_uplus, hya_uminus,hya_tposeBlk, hya_grv, hya_shk, hya_flat  
!!$  real, pointer, dimension(:,:,:,:) :: hy_flux,  hy_uplus, hy_uminus,hy_tposeBlk
!!$  real, pointer, dimension(:,:,:) :: hy_grv, hy_shk, hy_flat


  !! These are the variable that are offloaded
  integer :: hy_sizex,hy_sizey,hy_sizez
  integer :: hy_geometry
  integer, dimension(LOW:HIGH,MDIM) :: hy_lims

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
