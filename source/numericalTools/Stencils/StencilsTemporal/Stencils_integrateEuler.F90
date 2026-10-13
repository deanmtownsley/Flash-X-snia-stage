!!***
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
!!
!!***
#include "constants.h"

subroutine Stencils_integrateEulerScalar(phi,rhs,dt,lo,hi,iSource)
    implicit none
    real, dimension(:,:,:), intent(inout):: phi
    real, dimension(:,:,:), intent(in) :: rhs
    real, intent(in) :: dt
    integer,dimension(MDIM), intent(in) :: lo, hi
    real, intent(in) :: iSource

    phi(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) = &
           phi(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) + &
           dt*rhs(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) + &
           dt*iSource
    
    return
end subroutine Stencils_integrateEulerScalar

subroutine Stencils_integrateEulerArray(phi,rhs,dt,lo, hi,iSource)
    implicit none
    real, dimension(:,:,:), intent(inout):: phi
    real, dimension(:,:,:), intent(in) :: rhs
    real, intent(in) :: dt
    integer, dimension(MDIM), intent(in) :: lo, hi
    real, dimension(:,:,:), intent(in) :: iSource

    phi(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) = &
                   phi(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) + &
                   dt*rhs(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS)) + &
                   dt*iSource(lo(IAXIS):hi(IAXIS),lo(JAXIS):hi(JAXIS),lo(KAXIS):hi(KAXIS))
    return
end subroutine Stencils_integrateEulerArray
