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
subroutine Stencils_integrateAB2Scalar(phi,rhsNew,rhsOld,dt,ix1,ix2,jy1,jy2,kz1,kz2,iSource)
    implicit none
    real, dimension(:,:,:), intent(inout):: phi
    real, dimension(:,:,:), intent(in) :: rhsNew, rhsOld
    real, intent(in) :: dt
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: iSource

    phi(ix1:ix2,jy1:jy2,kz1:kz2) = phi(ix1:ix2,jy1:jy2,kz1:kz2) + &
                                   1.5*dt*rhsNew(ix1:ix2,jy1:jy2,kz1:kz2) - &
                                   0.5*dt*rhsOld(ix1:ix2,jy1:jy2,kz1:kz2) + &
                                        dt*iSource
    return
end subroutine Stencils_integrateAB2Scalar

subroutine Stencils_integrateAB2Array(phi,rhsNew,rhsOld,dt,ix1,ix2,jy1,jy2,kz1,kz2,iSource)
    implicit none
    real, dimension(:,:,:), intent(inout):: phi
    real, dimension(:,:,:), intent(in) :: rhsNew, rhsOld
    real, intent(in) :: dt
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, dimension(:,:,:), intent(in) :: iSource

    phi(ix1:ix2,jy1:jy2,kz1:kz2) = phi(ix1:ix2,jy1:jy2,kz1:kz2) + &
                                   1.5*dt*rhsNew(ix1:ix2,jy1:jy2,kz1:kz2) - &
                                   0.5*dt*rhsOld(ix1:ix2,jy1:jy2,kz1:kz2) + &
                                        dt*iSource(ix1:ix2,jy1:jy2,kz1:kz2)
    return
end subroutine Stencils_integrateAB2Array
