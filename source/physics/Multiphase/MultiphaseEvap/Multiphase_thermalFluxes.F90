!!***if* source/physics/Multiphase/MultiphaseEvap/Multiphase_thermalFluxes
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!!
!!***
!!REORDER(4): solnData

#include "constants.h"
#include "Multiphase.h"
#include "Simulation.h"

subroutine Multiphase_thermalFluxes()

  use Multiphase_data
  use Timers_interface,   ONLY : Timers_start, Timers_stop
  use Driver_interface,   ONLY : Driver_getNStep
  use Grid_interface,     ONLY : Grid_getTileIterator,Grid_releaseTileIterator,Grid_fillGuardCells
  use Grid_tile,          ONLY : Grid_tile_t
  use Grid_iterator,      ONLY : Grid_iterator_t
  use Stencils_interface, ONLY : Stencils_cnt_advectUpwind2d, Stencils_cnt_advectUpwind3d
  use mph_evapInterface,  ONLY : mph_phasedFluxes

!------------------------------------------------------------------------------------------------
  implicit none
  include "Flashx_mpi.h"
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  logical :: gcMask(NUNK_VARS+NDIM*NFACE_VARS)
  real, pointer, dimension(:,:,:,:) :: solnData
  integer :: ierr,i,j,k,ii
  real del(MDIM)
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
  integer TA(2),count_rate
  real*8  ET

!------------------------------------------------------------------------------------------------
  CALL SYSTEM_CLOCK(TA(1),count_rate)

  nullify(solnData)

  do ii=1,mph_extpIt
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%deltas(del)

#if NDIM < MDIM

     solnData(MFLX_VAR,:,:,:) = 0.0

     call Stencils_cnt_advectUpwind2d(solnData(MFLX_VAR,:,:,:),&
                                      solnData(HFLQ_VAR,:,:,:),&
                                      solnData(NRMX_VAR,:,:,:),&
                                      solnData(NRMY_VAR,:,:,:),&
                                      del(IAXIS),del(JAXIS),&
                                      GRID_ILO,GRID_IHI,&
                                      GRID_JLO,GRID_JHI)

     call mph_phasedFluxes(solnData(HFLQ_VAR,:,:,:),&
                           solnData(MFLX_VAR,:,:,:),&
                           solnData(DFUN_VAR,:,:,:),&
                           0.5*del(IAXIS),&
                           GRID_ILO,GRID_IHI,&
                           GRID_JLO,GRID_JHI,&
                           GRID_KLO,GRID_KHI)

     solnData(MFLX_VAR,:,:,:) = 0.0

      call Stencils_cnt_advectUpwind2d(solnData(MFLX_VAR,:,:,:),&
                                       solnData(HFGS_VAR,:,:,:),&
                                      -solnData(NRMX_VAR,:,:,:),&
                                      -solnData(NRMY_VAR,:,:,:),&
                                       del(IAXIS),del(JAXIS),&
                                       GRID_ILO,GRID_IHI,&
                                       GRID_JLO,GRID_JHI)

     call mph_phasedFluxes(solnData(HFGS_VAR,:,:,:),&
                           solnData(MFLX_VAR,:,:,:),&
                          -solnData(DFUN_VAR,:,:,:),&
                           0.5*del(IAXIS),&
                           GRID_ILO,GRID_IHI,&
                           GRID_JLO,GRID_JHI,&
                           GRID_KLO,GRID_KHI)

#else
 
     solnData(MFLX_VAR,:,:,:) = 0.0

     call Stencils_cnt_advectUpwind3d(solnData(MFLX_VAR,:,:,:),&
                                      solnData(HFLQ_VAR,:,:,:),&
                                      solnData(NRMX_VAR,:,:,:),&
                                      solnData(NRMY_VAR,:,:,:),&
                                      solnData(NRMZ_VAR,:,:,:),&
                                      del(IAXIS),del(JAXIS),del(KAXIS),&
                                      GRID_ILO,GRID_IHI,&
                                      GRID_JLO,GRID_JHI,&
                                      GRID_KLO,GRID_KHI)

     call mph_phasedFluxes(solnData(HFLQ_VAR,:,:,:),&
                           solnData(MFLX_VAR,:,:,:),&
                           solnData(DFUN_VAR,:,:,:),&
                           0.5*del(IAXIS),&
                           GRID_ILO,GRID_IHI,&
                           GRID_JLO,GRID_JHI,&
                           GRID_KLO,GRID_KHI)

     solnData(MFLX_VAR,:,:,:) = 0.0

      call Stencils_cnt_advectUpwind3d(solnData(MFLX_VAR,:,:,:),&
                                       solnData(HFGS_VAR,:,:,:),&
                                      -solnData(NRMX_VAR,:,:,:),&
                                      -solnData(NRMY_VAR,:,:,:),&
                                      -solnData(NRMZ_VAR,:,:,:),&
                                       del(IAXIS),del(JAXIS),del(KAXIS),&
                                       GRID_ILO,GRID_IHI,&
                                       GRID_JLO,GRID_JHI,&
                                       GRID_KLO,GRID_KHI)

     call mph_phasedFluxes(solnData(HFGS_VAR,:,:,:),&
                           solnData(MFLX_VAR,:,:,:),&
                          -solnData(DFUN_VAR,:,:,:),&
                           0.5*del(IAXIS),&
                           GRID_ILO,GRID_IHI,&
                           GRID_JLO,GRID_JHI,&
                           GRID_KLO,GRID_KHI)

     
#endif
    
      ! Release pointers:
      call tileDesc%releaseDataPtr(solnData, CENTER)
      call itor%next()
   end do
   call Grid_releaseTileIterator(itor)  

   gcMask = .FALSE.
   gcMask(HFLQ_VAR) = .TRUE.
   gcMask(HFGS_VAR) = .TRUE.
   call Grid_fillGuardCells(CENTER,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
 
  end do
  !
  !
  !
  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData,  CENTER)
     call tileDesc%deltas(del)

     solnData(MFLX_VAR,:,:,:) = (mph_Stefan*mph_invReynolds/mph_Prandtl)*&
                                (solnData(HFLQ_VAR,:,:,:)+mph_thcoGas*solnData(HFGS_VAR,:,:,:))  
 
     ! Release pointers:
     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
   end do
   call Grid_releaseTileIterator(itor)  

   gcMask = .FALSE.
   gcMask(MFLX_VAR) = .TRUE.
   call Grid_fillGuardCells(CENTER,ALLDIR,&
        maskSize=NUNK_VARS+NDIM*NFACE_VARS,mask=gcMask)
 
   CALL SYSTEM_CLOCK(TA(2),count_rate)
   ET=REAL(TA(2)-TA(1))/count_rate
   if (mph_meshMe .eq. MASTER_PE)  write(*,*) 'Multiphase thermalFluxes Time =',ET

end subroutine Multiphase_thermalFluxes
