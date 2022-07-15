!!****i** source/physics/Hydro/HydroMain/Spark/Hydro_funcs
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
!!  Hydro_funcs
!!
!! DESCRIPTION
!!
!!  Holds functions frequently accessed by Hydro.F90 for Spark hydrodynamic solver.
!! 
!!***
!!Reorder directive used by FLASH with --index-reorder flag at setup
!!Reorder(4): hy_starState,Uin, U, hy_fl[xyz], hy_fluxBuf[XYZ]

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

subroutine allocate_scr(blkLimits,blkLimitsGC)
  
  use Hydro_data, ONLY : hy_starState,  hy_fluxCorrect, hy_grav, hy_flx, hy_fly, hy_flz,&
       hy_tiny,hy_hybridRiemann,hy_C_hyp, &
       hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_cvisc, hy_del,hy_geometry, &
       hy_alphaGLM,hy_Vc,scratch_allocated
  use Hydro_data, ONLY : hydro_GPU_scratch, hy_uPlus, hy_uMinus,&
       hy_shck, hy_rope, hy_flux, hy_flat, hy_grv,hy_flat3d,hy_tmpState
  implicit none
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  integer :: max_edge, max_edge_y, max_edge_z
  
  max_edge = max(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 2,blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 2, &
       blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 2)
  max_edge_y = 1
  max_edge_z = 1
#if NDIM==2
  max_edge_y = max_edge
#elif NDIM==3
  max_edge_y = max_edge
  max_edge_z = max_edge
#endif
  ! print *, "max edge", max_edge
  !Construct arrays to hold fluxes used for solution update
  
  if (.NOT. allocated(hy_flx)) then
     allocate(hy_flx(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)+1,&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)+0*K2D,&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)+0*K3D))
     hy_flx = 0.
  endif
  
  if (.NOT. allocated(hy_fly)) then
     allocate(hy_fly(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)+0,&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)+1*K2D,&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)+0*K3D))
     hy_fly = 0.
  endif
  if (.NOT. allocated(hy_flz)) then
     allocate(hy_flz(NFLUXES,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS)+0,&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS)+0*K2D,&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)+1*K3D))
     hy_flz = 0.
  endif
  if (.NOT. allocated(hy_flat3d)) then
     allocate(hy_flat3d(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
     hy_flat3d = 0.
  endif
  if (.NOT. allocated(hy_Vc)) then
     allocate(hy_Vc(blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
     hy_Vc =0.
  end if
  
  
  
  !Gravity 
  if (.NOT. allocated(hy_grav)) then
     allocate(hy_grav(MDIM,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  endif
  
  if (.NOT. allocated(hy_starState)) then
     allocate(hy_starState(NUNK_VARS,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  endif
  
  if (.NOT. allocated(hy_tmpState)) then
     allocate(hy_tmpState(NUNK_VARS,blkLimitsGC(LOW,IAXIS):blkLimitsGC(HIGH,IAXIS),&
blkLimitsGC(LOW,JAXIS):blkLimitsGC(HIGH,JAXIS),&
blkLimitsGC(LOW,KAXIS):blkLimitsGC(HIGH,KAXIS)))
  endif

  !$omp target enter data map(alloc:hy_flat,hy_shck,hy_rope,hy_uMinus,hy_uPlus,hy_grv,hy_flux)

  call allocate_fxscr(blkLimits,blkLimitsGC)
  
end subroutine allocate_scr

subroutine deallocate_scr()
  use Hydro_data, ONLY : hy_starState,  hy_grav, hy_flx, hy_fly, hy_flz,&
       hy_Vc
  
  use Hydro_data, ONLY : &
       hy_flat3d,hy_tmpState

  if(allocated(hy_flx))deallocate(hy_flx)
  if(allocated(hy_fly))deallocate(hy_fly)
  if(allocated(hy_flz))deallocate(hy_flz)
  if(allocated(hy_flat3d))deallocate(hy_flat3d)
  if(allocated(hy_Vc))deallocate(hy_Vc)
  if(allocated(hy_grav))deallocate(hy_grav)
  if(allocated(hy_starState))deallocate(hy_starState)
  if(allocated(hy_tmpState))deallocate(hy_tmpState)
  call deallocate_fxscr()
end subroutine deallocate_scr

subroutine allocate_fxscr(blkLimits,blkLimitsGC)
  use Hydro_data, ONLY : hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ, &
       hy_eosData, hy_mfrac, hy_fluxCorrect
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  integer :: max_edge

  max_edge = max(blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS) + 2,blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS) + 2, &
       blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS) + 2)

  !Allocate size of flux buffers used for flux correction
  if (hy_fluxCorrect) then
     !allocate buffers here
     if (.NOT. allocated(hy_fluxBufX)) then 
        allocate(hy_fluxBufX(NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+1,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+0*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+0*K3D))
        hy_fluxBufX = 0.
     endif
     
     if (.NOT. allocated(hy_fluxBufY)) then 
        allocate(hy_fluxBufY(NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+0,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+1*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+0*K3D))
        hy_fluxBufY = 0.
     endif
     
     if (.NOT. allocated(hy_fluxBufZ)) then 
        allocate(hy_fluxBufZ(NFLUXES,blkLimits(LOW,IAXIS):blkLimits(HIGH,IAXIS)+0,&
blkLimits(LOW,JAXIS):blkLimits(HIGH,JAXIS)+0*K2D,&
blkLimits(LOW,KAXIS):blkLimits(HIGH,KAXIS)+1*K3D))
        hy_fluxBufZ = 0.
     endif
  endif
  !Set up one block's worth of local gravity.  Allocation allows for compatibility with Paramesh4 and AMRex

  if (.NOT. allocated(hy_eosData)) then
     allocate(hy_eosData(EOS_NUM*max_edge))
     hy_eosData = 0.
  endif
  
  if (.NOT. allocated(hy_mfrac)) then
     allocate(hy_mfrac(NSPECIES*max_edge))
     hy_mfrac = 0.
  endif
end subroutine allocate_fxscr


subroutine deallocate_fxscr()
  use Hydro_data, ONLY : hy_fluxBufX, hy_fluxBufY, hy_fluxBufZ, &
       hy_eosData, hy_mfrac, hy_fluxCorrect

  !Allocate size of flux buffers used for flux correction
  if (hy_fluxCorrect) then
     !allocate buffers here
     if(allocated(hy_fluxBufX))deallocate(hy_fluxBufX)
     if(allocated(hy_fluxBufY))deallocate(hy_fluxBufY)
     if(allocated(hy_fluxBufZ))deallocate(hy_fluxBufZ)
  end if

  if(allocated(hy_mfrac))deallocate(hy_mfrac)
  if(allocated(hy_eosData))deallocate(hy_eosData)
end subroutine deallocate_fxscr
