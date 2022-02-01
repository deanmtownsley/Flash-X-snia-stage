!!****f* source/Grid/Grid_pfftMapToInput
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
!!   Grid_pfftMapToInput
!!
!! SYNOPSIS
!!
!!   Grid_pfftMapToInput(integer(IN) :: gridVar,
!!                 real(INOUT) :: pfft_inArray)
!!
!! DESCRIPTION 
!!
!!  Takes the data correnspoding to 
!!  the variable gridVar in the mesh data structure 
!!  and redistributes to make it campatible with Pfft requirements
!!  based upon thh map determined by the routine Grid_pfftInit.
!!
!! ARGUMENTS
!!
!!  gridVar          - variable on the mesh on which pfft is to be applies
!!  pfft_inArray     - array that is input to pfft
!!
!! NOTE 
!!
!!  Users must make sure that Grid_pfftInput has been called correctly
!!  before calling this routine.
!!
!!***
subroutine Grid_pfftMapToInput(gridVar,pfft_inArray)
  implicit none
  integer,intent(IN) :: gridVar
  real,dimension(:),intent(inout) :: pfft_inArray
  pfft_inArray=0.0
end subroutine Grid_pfftMapToInput
