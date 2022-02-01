!!
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
!!**
subroutine Stencils_advectCentral2d(rhs,phi,u,v,dx,dy,ix1,ix2,jy1,jy2,center,facex,facey)

  implicit none

  !----Imported variables
  real, dimension(:,:,:), intent(inout):: rhs
  real, dimension(:,:,:), intent(in) :: phi,u,v
  real, intent(in) :: dx,dy
  integer, intent(in) :: ix1,ix2,jy1,jy2
  integer, intent(in) :: center,facex,facey

  !----Local variables-----
  integer:: i, j
  real:: dx1, dy1
  real:: uplus,umins,vplus,vmins
  real:: phi_uplus,phi_umins,phi_vplus,phi_vmins
  integer, parameter :: kz1 = 1

  ! grid spacings
  dx1 = 1.0/dx
  dy1 = 1.0/dy

  do j = jy1,jy2
     do i = ix1,ix2

        umins = center*u(i,j,kz1) + &
                facex*0.5*(u(i,j,kz1) + u(i-1,j,kz1)) + &
                facey*0.5*(u(i,j,kz1) + u(i,j-1,kz1))

        uplus = center*u(i+1,j,kz1) + &
                facex*0.5*(u(i+1,j,kz1)+u(i,j,kz1)) + &
                facey*0.5*(u(i+1,j,kz1)+u(i+1,j-1,kz1))

        vmins = center*v(i,j,kz1) + &
                facex*0.5*(v(i,j,kz1) + v(i-1,j,kz1)) + &
                facey*0.5*(v(i,j,kz1) + v(i,j-1,kz1))

        vplus = center*v(i,j+1,kz1) + &
                facex*0.5*(v(i,j+1,kz1) + v(i-1,j+1,kz1)) + &
                facey*0.5*(v(i,j+1,kz1) + v(i,j,kz1))

        phi_uplus = (phi(i+1,j,kz1) + phi(i,j,kz1))*0.5
        phi_umins = (phi(i,j,kz1)   + phi(i-1,j,kz1))*0.5
        phi_vplus = (phi(i,j+1,kz1) + phi(i,j,kz1))*0.5
        phi_vmins = (phi(i,j,kz1)   + phi(i,j-1,kz1))*0.5

        ! calculate RHS for u-momentum
        rhs(i,j,kz1) =    rhs(i,j,kz1)                                  &
                          - (phi_uplus*uplus - phi_umins*umins)*dx1     &! advection term
                          - (phi_vplus*vplus - phi_vmins*vmins)*dy1                          

      enddo
   enddo
   return
end subroutine Stencils_advectCentral2d
