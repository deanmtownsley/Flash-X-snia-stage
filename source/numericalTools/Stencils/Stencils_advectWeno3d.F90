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
subroutine Stencils_advectWeno3d(rhs,phi,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,&
                                  center,facex,facey,facez)     
  implicit none
  real, dimension(:,:,:), intent(inout):: rhs
  real, dimension(:,:,:), intent(in) :: phi,u,v,w
  real, intent(in) :: dx,dy,dz
  integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
  integer, intent(in) :: center,facex,facey,facez
end subroutine Stencils_advectWeno3d
