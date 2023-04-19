!!****h* source/physics/Eos/EosMain/WeakLib/Eos_wlInterface
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!! This is the header file for the WeakLib unit that defines its
!! public interfaces.
!!***


Module Eos_wlInterface

  implicit none

  interface
    subroutine Eos_wlOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xdedt,xCs2,xXp,xXn,xXa,xXh,xAbar,xVar,varID,mode)
    implicit none
    real, intent(INOUT) :: xDens, xYe
    real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
    integer, intent(IN) :: mode, varID
    real, intent(OUT) :: xXp, xXn, xXa,xXh,xdedt,xCs2,xVar,xAbar
    real :: xZbar,xMu_e,xMu_n,xMu_p,xMuhat
    end subroutine Eos_wlOneZone
  end interface

  interface
     subroutine Eos_wlDetectBounce (postBounce,bounceTime,centralDens,centralEntr)
       implicit none
       logical, intent(OUT) :: postBounce
       real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr
     end subroutine Eos_wlDetectBounce
  end interface

end Module Eos_wlInterface
  
