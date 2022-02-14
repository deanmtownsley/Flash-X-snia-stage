!!****if* /nfs2/kweide/projects/FLASH5/source/Grid/GridMain/AMR/Paramesh4/gr_wdInit
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
!!  gr_wdInit
!!
!! SYNOPSIS
!!
!!  call gr_wdInit()
!!
!! DESCRIPTION
!!
!!  Initializations for work distribution.
!!  
!!  Initialize runtime parameters that control the load distribution
!!  for PARAMESH.
!!
!! ARGUMENTS
!!
!!   No arguments
!!
!! NOTES
!!
!!   Runtime parameters should probably be named gr_wd* instead of gr_bt*.
!!
!!***

subroutine gr_wdInit()
  use RuntimeParameters_interface, ONLY: RuntimeParameters_get
  use Grid_data, ONLY: gr_meshMe
  use tree,      ONLY: gr_btSortByWork,       gr_btCustomWork,     &
                       gr_btWorkDefaultLeaf,  gr_btWorkDefaultPar, &
                       gr_btWorkBoundsLeaf,   gr_btWorkBoundsPar,  &
                       gr_btWorkChildScaling,                      &
                       gr_btExchangeWork

#include "constants.h"

  implicit none

! Set runtime parameters controlling load distribution
  call RuntimeParameters_get("gr_btSortByWork",gr_btSortByWork)
  call RuntimeParameters_get("gr_btCustomWork",gr_btCustomWork)
  call RuntimeParameters_get("gr_btWorkDefaultLeaf",gr_btWorkDefaultLeaf)
  call RuntimeParameters_get("gr_btWorkDefaultPar",gr_btWorkDefaultPar)
  call RuntimeParameters_get("gr_btWorkLBLeaf",gr_btWorkBoundsLeaf(LOW))
  call RuntimeParameters_get("gr_btWorkLBPar",gr_btWorkBoundsPar(LOW))
  call RuntimeParameters_get("gr_btWorkUBLeaf",gr_btWorkBoundsLeaf(HIGH))
  call RuntimeParameters_get("gr_btWorkUBPar",gr_btWorkBoundsPar(HIGH))
  call RuntimeParameters_get("gr_btWorkChildScaling",gr_btWorkChildScaling)
  call RuntimeParameters_get("gr_btExchangeWork",gr_btExchangeWork)
  if(gr_btCustomWork.AND..NOT.gr_btSortByWork) then
     if (gr_meshMe == MASTER_PE) &
          print *,"Warning: gr_btCustomWork &
                  &only has effect if gr_btSortByWork is True"
  end if
end subroutine gr_wdInit
