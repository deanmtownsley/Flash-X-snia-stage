!!****if* source/Grid/GridMain/AMR/Paramesh4/Grid_finalize
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
!!
!!  Grid_finalize
!!
!!
!! SYNOPSIS
!!
!!  call Grid_finalize()
!!
!!
!! DESCRIPTION
!!
!!  Deallocates memory allocated in the Grid Unit to prepare for shutdowns
!!
!!***

#include "Simulation.h"

subroutine Grid_finalize()

  use Grid_data,      ONLY : gr_delta
  use gr_bcInterface, ONLY : gr_bcFinalize
  use gr_ptInterface, ONLY : gr_ptFinalize
  use gr_specificData, ONLY : gr_nToLeft, gr_gid, gr_gsurr_blks, &
       scratch,scratch_ctr,scratch_facevarx,scratch_facevary,scratch_facevarz, &
       gr_flxx,gr_flxy,gr_flxz, gr_xflx,gr_yflx,gr_zflx
#ifdef FLASH_HYDRO_UNSPLIT
# if N_DIM >= 2
  use gr_specificData, ONLY : gr_xflx_yface,gr_yflx_xface
#  if N_DIM == 3
  use gr_specificData, ONLY : gr_xflx_zface,gr_yflx_zface, &
                              gr_zflx_xface,gr_zflx_yface
#  endif
# endif
#endif

  use gr_sbInterface, ONLY: gr_sbFinalize
  implicit none

  if(allocated(gr_gid))deallocate(gr_gid)
  if(allocated(gr_nToLeft))deallocate(gr_nToLeft)
#ifdef FLASH_GRID_PARAMESH3OR4
  if(allocated(gr_gsurr_blks))deallocate(gr_gsurr_blks)
#endif

  call gr_ptFinalize()
  call gr_solversFinalize()
  call gr_bcFinalize()
!  call gr_sbFinalize()

#ifndef BSS_GRID_ARRAYS
  ! Cf. gr_initSpecific
  deallocate(scratch)
  deallocate(scratch_ctr)
  deallocate(scratch_facevarx)
  deallocate(scratch_facevary)
  deallocate(scratch_facevarz)
  deallocate(gr_flxx)
  deallocate(gr_flxy)
  deallocate(gr_flxz)
  deallocate(gr_xflx)
  deallocate(gr_yflx)
  deallocate(gr_zflx)
# ifdef FLASH_HYDRO_UNSPLIT
#  if N_DIM >= 2
  deallocate(gr_xflx_yface)
  deallocate(gr_yflx_xface)
#   if N_DIM == 3
  deallocate(gr_xflx_zface)
  deallocate(gr_yflx_zface)
  deallocate(gr_zflx_xface)
  deallocate(gr_zflx_yface)
#   endif
#  endif
# endif
#endif
  deallocate(gr_delta)

  call Paramesh_finalize()

end subroutine Grid_finalize
