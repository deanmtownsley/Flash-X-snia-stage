!!****if* source/Grid/GridMain/AMR/Amrex/gr_initSpecific
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
!!  gr_initSpecific
!!
!! SYNOPSIS
!!
!!  call gr_initSpecific()
!!
!!
!! DESCRIPTION
!!  Initialize some implementation-specific Grid data
!!
!!  This routine should initialize data in the gr_specificData module.
!!
!! ARGUMENTS
!!  none
!!
!!***


subroutine gr_initSpecific()
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Grid_data, ONLY : gr_meshNumProcs
  use gr_specificData, ONLY : gr_nToLeft, &
                              gr_bndGCFillNeedsPrimitiveVars

  implicit none


  !Initialize grid arrays used by IO
  allocate(gr_nToLeft(0:gr_meshNumProcs-1))

  call RuntimeParameters_get("gr_bndGCFillNeedsPrimitiveVars", gr_bndGCFillNeedsPrimitiveVars)

end subroutine gr_initSpecific

