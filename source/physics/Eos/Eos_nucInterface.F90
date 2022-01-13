Module Eos_nucInterface
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.

  implicit none

  interface
     subroutine Eos_nucOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xdedt,xCs2,xXp,xXn,xXa,xXh, xVar,varID,mode)
       implicit none
       real, intent(INOUT) :: xDens
       real, intent(IN)    :: xYe
       real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
       real, intent(OUT)   :: xVar
       integer, intent(IN) :: mode, varID
       real, intent(OUT) :: xXp, xXn, xXa,xXh,xdedt,xCs2
     end subroutine Eos_nucOneZone
  end interface

  interface
     subroutine Eos_nucDetectBounce (postBounce,bounceTime,centralDens,centralEntr)
       implicit none
       logical, intent(OUT) :: postBounce
       real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr
     end subroutine Eos_nucDetectBounce
  end interface

end Module Eos_nucInterface
  
