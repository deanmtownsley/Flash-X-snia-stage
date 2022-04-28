!!****if* source/Grid/GridBoundaryConditions/Flash3HSE/gr_bcHseData
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
!!  gr_bcHseData
!!
!! SYNOPSIS
!!
!!  use gr_bcHseData
!!
!!
!!***

Module gr_bcHseData
  integer,parameter :: HSE_FORWARD = 1
  integer,parameter :: HSE_BACKWARD = 2
  integer,parameter :: HSE_CONSTENTR = 3
  integer,parameter :: HSE_CONSTTEMP = 4
  integer,parameter :: HSE_SETTEMP = 5

  character(len=80), save :: gr_bcHseGravDirec
  integer, save :: gr_bcHseDirection

  real, save :: gr_bcHseGravConst

end Module gr_bcHseData
