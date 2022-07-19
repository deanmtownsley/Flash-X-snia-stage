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
!!Reorder(4): hya_starState,Uin, U, hy_fl[xyz], hy_fluxBuf[XYZ]

#include "Simulation.h"
#include "constants.h"
#include "Spark.h"
#include "Eos.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

subroutine allocate_scr(blkLimits,blkLimitsGC)
  
  use Hydro_data, ONLY :  hy_fluxCorrect, hya_grav, hya_flx, hya_fly, hya_flz,&
       hy_tiny,hy_hybridRiemann,hy_C_hyp, &
       hy_smalldens, hy_smallE, hy_smallpres, hy_smallX, hy_cvisc, hy_del,hy_geometry, &
       hy_alphaGLM

  use Hydro_data, ONLY : hya_starState,  hya_uPlus, hya_uMinus,hya_Vc,&
       hya_shck, hya_rope, hya_flux, hya_flat, hya_grv,hya_flat3d,hya_tmpState

  implicit none
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  
  integer :: max_edge, max_edge_y, max_edge_z,space
  
  max_edge = MEDGE+2
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
  space=max_edge*max_edge_y*max_edge_z
  if (.NOT. allocated(hya_flux)) then
     allocate(hya_flux(NFLUXES*space))
     hya_flux=0.
  end if

  if (.NOT. allocated(hya_rope)) then
     allocate(hya_rope((1+NRECON)*space))
     hya_rope=0.
  end if

  if (.NOT. allocated(hya_uPlus)) then
     allocate(hya_uPlus((1+NRECON)*space))
     hya_uPlus=0.
  end if

  if (.NOT. allocated(hya_uMinus)) then
     allocate(hya_uMinus((1+NRECON)*space))
     hya_uMinus=0.
  end if

  if (.NOT. allocated(hya_flat)) then
     allocate(hya_flat(space))
     hya_flat=0.
  end if
  
  if (.NOT. allocated(hya_shck)) then
     allocate(hya_shck(space))
     hya_shck=0.
  end if
  
  if (.NOT. allocated(hya_grv)) then
     allocate(hya_grv(space))
     hya_grv=0.
  end if

  if (.NOT. allocated(hya_starState)) then
     allocate(hya_starState(NUNK_VARS*space))
  endif
  
  if (.NOT. allocated(hya_tmpState)) then
     allocate(hya_tmpState(NUNK_VARS*space))
  endif

  if (.NOT. allocated(hya_flat3d)) then
     allocate(hya_flat3d(space))
     hya_flat3d = 0.
  endif
  if (.NOT. allocated(hya_Vc)) then
     allocate(hya_Vc(space))
     hya_Vc =0.
  end if
  
    
  !Gravity 
  if (.NOT. allocated(hya_grav)) then
     allocate(hya_grav(MDIM*space))
  endif

  
  
  if (.NOT. allocated(hya_flx)) then
     allocate(hya_flx(NFLUXES*space))
     hya_flx = 0.
  endif
  
  if (.NOT. allocated(hya_fly)) then
     allocate(hya_fly(NFLUXES*space))
     hya_fly = 0.
  endif
  if (.NOT. allocated(hya_flz)) then
     allocate(hya_flz(NFLUXES*Space))
     hya_flz = 0.
  endif
  
  
  call allocate_fxscr(blkLimits,blkLimitsGC)
  
end subroutine allocate_scr

subroutine deallocate_scr()
  use Hydro_data, ONLY : hya_starState,  hya_grav, hya_flx, hya_fly, hya_flz,&
       hya_Vc,  hya_flat3d,hya_tmpState
  

  if(allocated(hya_flx))deallocate(hya_flx)
  if(allocated(hya_fly))deallocate(hya_fly)
  if(allocated(hya_flz))deallocate(hya_flz)
  if(allocated(hya_flat3d))deallocate(hya_flat3d)
  if(allocated(hya_Vc))deallocate(hya_Vc)
  if(allocated(hya_grav))deallocate(hya_grav)
  if(allocated(hya_starState))deallocate(hya_starState)
  if(allocated(hya_tmpState))deallocate(hya_tmpState)
  call deallocate_fxscr()
end subroutine deallocate_scr

subroutine allocate_fxscr(blkLimits,blkLimitsGC)
  use Hydro_data, ONLY : hya_fluxBufX, hya_fluxBufY, hya_fluxBufZ, &
       hy_eosData, hy_mfrac, hy_fluxCorrect
  integer,dimension(LOW:HIGH,MDIM),intent(IN) :: blkLimits, blkLimitsGC
  integer :: max_edge, max_edge_y, max_edge_z,space
  
  max_edge = MEDGE+2
  max_edge_y = 1
  max_edge_z = 1
#if NDIM==2
  max_edge_y = max_edge
#elif NDIM==3
  max_edge_y = max_edge
  max_edge_z = max_edge
#endif

  space=max_edge*max_edge_y*max_edge_z

  !Allocate size of flux buffers used for flux correction
  if (hy_fluxCorrect) then
     !allocate buffers here
     if (.NOT. allocated(hya_fluxBufX)) then 
        allocate(hya_fluxBufX(NFLUXES*space))
        hya_fluxBufX = 0.
     endif
     
     if (.NOT. allocated(hya_fluxBufY)) then 
        allocate(hya_fluxBufY(NFLUXES*space))
        hya_fluxBufY = 0.
     endif
     
     if (.NOT. allocated(hya_fluxBufZ)) then 
        allocate(hya_fluxBufZ(NFLUXES*space))
        hya_fluxBufZ = 0.
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
  use Hydro_data, ONLY : hya_fluxBufX, hya_fluxBufY, hya_fluxBufZ, &
       hy_eosData, hy_mfrac, hy_fluxCorrect

  !Allocate size of flux buffers used for flux correction
  if (hy_fluxCorrect) then
     !allocate buffers here
     if(allocated(hya_fluxBufX))deallocate(hya_fluxBufX)
     if(allocated(hya_fluxBufY))deallocate(hya_fluxBufY)
     if(allocated(hya_fluxBufZ))deallocate(hya_fluxBufZ)
  end if

  if(allocated(hy_mfrac))deallocate(hy_mfrac)
  if(allocated(hy_eosData))deallocate(hy_eosData)
end subroutine deallocate_fxscr
