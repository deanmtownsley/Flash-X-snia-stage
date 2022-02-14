subroutine ins_omgToCenter(uf,vf,wf,omgm,ix1,ix2,jy1,jy2,kz1,kz2,dx,dy,dz)
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

#include "constants.h"
#include "Simulation.h"

    implicit none

    !-------Arguments---------!
    real, dimension(:,:,:), intent(in)  :: uf,vf,wf
    real, dimension(:,:,:), intent(out) :: omgm
    integer, intent(in) :: ix1,jy1,kz1
    integer, intent(in) :: ix2,jy2,kz2
    real, intent(in) :: dx,dy,dz

    !-------Local Variables-----!
    integer :: i,j,k
    real, dimension(ix2,jy2,kz2) :: omgvertz,omgz
#if NDIM == MDIM
    real, dimension(ix2,jy2,kz2) :: omgvertx,omgx
    real, dimension(ix2,jy2,kz2) :: omgverty,omgy
#endif

    ! Compute omgvert on z direction, minimum number of NGUARD == 2.
    omgvertz(ix1+1:ix2,jy1+1:jy2,kz1:kz2) = (vf(ix1+1:ix2,jy1+1:jy2,kz1:kz2)-vf(ix1:ix2-1,jy1+1:jy2,kz1:kz2))/dx + &
                                            (uf(ix1+1:ix2,jy1+1:jy2,kz1:kz2)-uf(ix1+1:ix2,jy1:jy2-1,kz1:kz2))/dy

      
    omgz(ix1+1:ix2-1,jy1+1:jy2-1,kz1:kz2) = &
    0.25*(omgvertz(ix1+1:ix2-1,jy1+1:jy2-1,kz1:kz2)+ &
          omgvertz(ix1+2:ix2  ,jy1+1:jy2-1,kz1:kz2)+ &
          omgvertz(ix1+1:ix2-1,jy1+2:jy2  ,kz1:kz2)+ &
          omgvertz(ix1+2:ix2  ,jy1+2:jy2  ,kz1:kz2))

   omgm(ix1+1:ix2-1,jy1+1:jy2-1,kz1:kz2) = omgz(ix1+1:ix2-1,jy1+1:jy2-1,kz1:kz2)

end subroutine ins_omgToCenter
