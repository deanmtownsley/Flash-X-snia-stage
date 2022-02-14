!!***if* source/physics/HeatAD/HeatADAdvection/HeatAD_advection
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
!!
!!
!!***
!!REORDER(4): solnData,face[xyz]Data

#include "constants.h"
#include "HeatAD.h"
#include "Simulation.h"   

subroutine HeatAD_advection()

   use HeatAD_data
   use Timers_interface,    ONLY : Timers_start, Timers_stop
   use Driver_interface,    ONLY : Driver_getNStep
   use Grid_interface,      ONLY : Grid_getTileIterator,Grid_releaseTileIterator
   use Grid_tile,           ONLY : Grid_tile_t
   use Grid_iterator,       ONLY : Grid_iterator_t
   use Stencils_interface,  ONLY : Stencils_advectWeno2d,Stencils_advectWeno3d

!--------------------------------------------------------------------------------------------
   implicit none
   include"Flashx_mpi.h"
   real ::  del(MDIM)
   integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
   real, pointer, dimension(:,:,:,:) :: solnData,facexData,faceyData,facezData
   integer TA(2),count_rate,ierr
   real*8  ET
   type(Grid_tile_t) :: tileDesc
   type(Grid_iterator_t) :: itor

!---------------------------------------------------------------------------------------------
   CALL SYSTEM_CLOCK(TA(1),count_rate)

   nullify(solnData,facexData,faceyData,facezData)

   call Grid_getTileIterator(itor, nodetype=LEAF)
   do while(itor%isValid())
       call itor%currentTile(tileDesc)
       call tileDesc%getDataPtr(solnData,  CENTER)
       call tileDesc%deltas(del)

#if NDIM == MDIM
       call tileDesc%getDataPtr(facexData,  FACEX)
       call tileDesc%getDataPtr(faceyData,  FACEY)
       call tileDesc%getDataPtr(facezData,  FACEZ)

       call Stencils_advectWeno3d(solnData(RHST_VAR,:,:,:),&
                                  solnData(TEMP_VAR,:,:,:),&
                                  facexData(ht_iVelFVar,:,:,:),&
                                  faceyData(ht_iVelFVar,:,:,:),&
                                  facezData(ht_iVelFvar,:,:,:),&
                                  del(DIR_X),del(DIR_Y),del(DIR_Z),&
                                  GRID_ILO, GRID_IHI, &
                                  GRID_JLO, GRID_JHI, &
                                  GRID_KLO, GRID_KHI, &
                                  center=.true.,facex=.false.,facey=.false.,facez=.false.)

       call tileDesc%releaseDataPtr(facexData, FACEX)
       call tileDesc%releaseDataPtr(faceyData, FACEY)
       call tileDesc%releaseDataPtr(facezData, FACEZ)

#else

       call tileDesc%getDataPtr(facexData,  FACEX)
       call tileDesc%getDataPtr(faceyData,  FACEY)

       call Stencils_advectWeno2d(solnData(RHST_VAR,:,:,:),&
                                  solnData(TEMP_VAR,:,:,:),&
                                  facexData(ht_iVelFVar,:,:,:),&
                                  faceyData(ht_iVelFVar,:,:,:),&
                                  del(DIR_X),del(DIR_Y),&
                                  GRID_ILO, GRID_IHI, &
                                  GRID_JLO, GRID_JHI, &
                                  center=.true.,facex=.false.,facey=.false.)

       call tileDesc%releaseDataPtr(facexData, FACEX)
       call tileDesc%releaseDataPtr(faceyData, FACEY)

#endif
       call tileDesc%releaseDataPtr(solnData, CENTER)
       call itor%next()
  end do
  call Grid_releaseTileIterator(itor)  

  CALL SYSTEM_CLOCK(TA(2),count_rate)
  ET=REAL(TA(2)-TA(1))/count_rate
  if (ht_meshMe .eq. MASTER_PE)  write(*,*) 'Total Heat AD Advection Time =',ET

end subroutine HeatAD_advection
