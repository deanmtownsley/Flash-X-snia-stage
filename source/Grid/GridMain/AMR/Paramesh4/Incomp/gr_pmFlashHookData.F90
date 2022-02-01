!!****ih* source/Grid/GridMain/paramesh/Incomp/gr_pmFlashHookData
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
!!  gr_pmFlashHookData
!!
!! SYNOPSIS
!!
!!  use gr_pmFlashHookData
!!
!! DESCRIPTION
!!
!!  Module for some data items used by FLASH to customize PARAMESH.
!!
!!***

! Modification history:
!     Created   October 2016  KW

module gr_pmFlashHookData

  implicit none
      
  logical,parameter :: gr_pmAlwaysFillFcGcAtDomainBC = .TRUE.

end module gr_pmFlashHookData
