!!****if* source/Grid/localAPI/gr_sbCreateGroups
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
!!  gr_sbCreateGroups
!!
!! SYNOPSIS
!!
!!  gr_sbCreateGroups()
!!  
!! DESCRIPTION 
!!  
!!  This routine is called from Grid_initDomain. It creates communicators
!!  for each solid body and identifies the processor that overlaps the most
!!  with the solid body as the master processor.
!!
!!
!! ARGUMENTS 
!!
!!***


Subroutine gr_sbCreateGroups()
  implicit none
End Subroutine gr_sbCreateGroups
