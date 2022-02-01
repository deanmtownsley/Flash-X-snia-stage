!!****if* source/Particles/localAPI/pt_findTagOffset
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
!!    pt_findTagOffset
!!
!! SYNOPSIS
!!
!!    pt_findTagOffset()
!!
!! DESCRIPTION
!!    Add a unique tag each to particles generated during evolution.
!!    The algorithm first finds out the sum of number of particles 
!!    being added to all the processors to the left of MyPE. It then
!!    adds pt_startTagNumber (the largest tag number in use) to it.
!!    This number acts as the offset for the tag numbers being assigned
!!    to the newly added particles.
!!
!! NOTES
!!    This method of tag generation will work for up to 10^14 
!!    particles in a simulation
!!
!!
!!
!!
!!***

!!#define DEBUG_PARTICLES

subroutine pt_findTagOffset(newCount,tagOffset)

  implicit none

  integer, intent(IN) :: newCount
  integer, intent(OUT) :: tagOffset
  
  tagOffset=0

  return
end subroutine pt_findTagOffset
