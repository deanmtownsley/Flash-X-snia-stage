!!****f* source/physics/Eos/Eos_nucDetectBounce
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!! NAME
!!  
!!  Eos_nucDetectBounce 
!!
!!
!! SYNOPSIS
!! 
!!  call Eos_nucDetectBounce(logical(OUT)  :: postbounce,
!!                           real, optional(OUT)  :: bouncetime,
!!                           real, optional(OUT)  :: centraldens,
!!                           real, optional(OUT)  :: centralentr)
!!
!!  
!!  
!! DESCRIPTION
!!  This routine determines if collapse has proceeded to the point of 
!!  core bounce, as determined by the maximum density.
!!
!! ARGUMENTS
!!
!!   postbounce : 
!!
!!   bouncetime : 
!!
!!   centraldens : 
!!
!!   centralentr : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine Eos_nucDetectBounce(postBounce,bounceTime,centralDens,centralEntr)

  implicit none
  logical, intent(OUT) :: postBounce
  real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr

  postBounce = .FALSE.
  if (present(bounceTime)) bounceTime = 0.0
  if (present(centralDens)) centralDens = 0.0
  if (present(centralEntr)) centralEntr = 0.0

  return
end subroutine Eos_nucDetectBounce
