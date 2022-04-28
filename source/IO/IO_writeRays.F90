!!****f* source/IO/IO_writeRays
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
!!  IO_writeRays
!!
!! SYNOPSIS
!!
!!  call IO_writeRays(integer(in) :: numrays,
!!                    integer(in) :: raytags,
!!                    real(in) :: posbuffer,
!!                    real(in) :: powerbuffer,
!!                    integer(in) :: numpos)
!!
!! DESCRIPTION
!! 
!! Write rays 
!!
!! ARGUMENTS
!!
!!   numrays : number of rays  
!!
!!   raytags : tags of rays 
!!
!!   posbuffer : position buffer
!!
!!   powerbuffer : power buffer
!!
!!   numpos : 
!!
!!
!!
!!***

subroutine IO_writeRays(numRays, rayTags, posBuffer, powerBuffer, numPos)
  implicit none

  integer, intent(in) :: numRays
  integer, intent(in) :: rayTags(:)
  real, intent(in) :: posBuffer(:,:,:)
  real, intent(in) :: powerBuffer(:,:)
  integer, intent(in) :: numPos(:)

end subroutine IO_writeRays
