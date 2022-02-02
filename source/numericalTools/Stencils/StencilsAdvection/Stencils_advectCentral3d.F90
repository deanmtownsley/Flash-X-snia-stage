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
subroutine Stencils_advectCentral3d(rhs,phi,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,&
                                     center,facex,facey,facez)
  implicit none

  !----Imported variables
  real, dimension(:,:,:), intent(inout):: rhs
  real, dimension(:,:,:), intent(in) :: phi,u,v,w
  real, intent(in) :: dx,dy,dz
  integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
  integer, intent(in) :: center,facex,facey,facez

  integer:: i, j, k
  real:: dx1, dy1, dz1
  real:: uplus, umins, vplus, vmins, wplus, wmins
  real:: phi_uplus, phi_umins, phi_vplus, phi_vmins, phi_wplus, phi_wmins

  ! grid spacings
  dx1 = 1.0/dx
  dy1 = 1.0/dy
  dz1 = 1.0/dz

  do k = kz1,kz2
     do j = jy1,jy2
        do i = ix1,ix2

           umins = center*u(i,j,k) + &
                   facex*0.5*(u(i,j,k) + u(i-1,j,k)) + &
                   facey*0.5*(u(i,j,k) + u(i,j-1,k)) + &
                   facez*0.5*(u(i,j,k) + u(i,j,k-1))

           uplus = center*u(i+1,j,k) + &
                   facex*0.5*(u(i+1,j,k)+u(i,j,k)) + &
                   facey*0.5*(u(i+1,j,k)+u(i+1,j-1,k)) + &
                   facez*0.5*(u(i+1,j,k)+u(i+1,j,k-1))

           vmins = center*v(i,j,k) + &
                   facex*0.5*(v(i,j,k) + v(i-1,j,k)) + &
                   facey*0.5*(v(i,j,k) + v(i,j-1,k)) + &
                   facez*0.5*(v(i,j,k) + v(i,j,k-1))

           vplus = center*v(i,j+1,k) + &
                   facex*0.5*(v(i,j+1,k) + v(i-1,j+1,k)) + &
                   facey*0.5*(v(i,j+1,k) + v(i,j,k)) + &
                   facez*0.5*(v(i,j+1,k) + v(i,j+1,k-1))

           wmins = center*w(i,j,k) + &
                   facex*0.5*(w(i,j,k) + w(i-1,j,k)) + &
                   facey*0.5*(w(i,j,k) + w(i,j-1,k)) + &
                   facez*0.5*(w(i,j,k) + w(i,j,k-1))

           wplus = center*w(i,j,k+1) + &
                   facex*0.5*(w(i,j,k+1) + w(i-1,j,k+1)) + &
                   facey*0.5*(w(i,j,k+1) + w(i,j-1,k+1)) + &
                   facez*0.5*(w(i,j,k+1) + w(i,j,k))

           phi_uplus  = (phi(i+1,j  ,k  ) + phi(i  ,j  ,k  ))*0.5
           phi_umins  = (phi(i  ,j  ,k  ) + phi(i-1,j  ,k  ))*0.5

           phi_vplus  = (phi(i  ,j+1,k  ) + phi(i  ,j  ,k  ))*0.5
           phi_vmins  = (phi(i  ,j  ,k  ) + phi(i  ,j-1,k  ))*0.5

           phi_wplus  = (phi(i  ,j  ,k+1) + phi(i  ,j  ,k  ))*0.5
           phi_wmins  = (phi(i  ,j  ,k  ) + phi(i  ,j  ,k-1))*0.5

           ! calculate RHS for u-momentum
           rhs(i,j,k)  = rhs(i,j,k)                               &                              
                       - (uplus*phi_uplus - umins*phi_umins)*dx1  &! advection term
                       - (vplus*phi_vplus - vmins*phi_vmins)*dy1  &
                       - (wplus*phi_wplus - wmins*phi_wmins)*dz1
        enddo
     enddo
  enddo
  return
end subroutine Stencils_advectCentral3d
