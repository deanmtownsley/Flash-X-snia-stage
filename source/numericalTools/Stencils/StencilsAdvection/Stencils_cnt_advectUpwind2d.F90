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
subroutine Stencils_cnt_advectUpwind2d(rhs,phi,u,v,dx,dy,ix1,ix2,jy1,jy2)

  implicit none

  !----Imported variables
  real, dimension(:,:,:), intent(inout):: rhs
  real, dimension(:,:,:), intent(in) :: phi,u,v
  real, intent(in) :: dx,dy
  integer, intent(in) :: ix1,ix2,jy1,jy2

  !----Local variables
  integer :: i,j,k
  real :: u_plus, u_mins, v_plus, v_mins, u_conv, v_conv
  real :: phix_plus, phix_mins, phiy_plus, phiy_mins

  k = 1

  do j=jy1,jy2
     do i=ix1,ix2

     u_conv = u(i,j,k)
     v_conv = v(i,j,k)

     u_plus = max(u_conv, 0.)
     u_mins = min(u_conv, 0.)

     v_plus = max(v_conv, 0.)
     v_mins = min(v_conv, 0.)

     phix_plus = phi(i+1,j,k)-phi(i,j,k)
     phix_mins = phi(i,j,k)-phi(i-1,j,k)

     phiy_plus = phi(i,j+1,k)-phi(i,j,k)
     phiy_mins = phi(i,j,k)-phi(i,j-1,k)

     rhs(i,j,k) = rhs(i,j,k) &
                             - (1./dx) * (u_plus*phix_mins + u_mins*phix_plus)&
                             - (1./dy) * (v_plus*phiy_mins + v_mins*phiy_plus)


    end do
  end do 

  return

end subroutine Stencils_cnt_advectUpwind2d
