!!***if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/outflow/sim_outflowVelBlk
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!!***

#include "constants.h"
#include "Simulation.h"

subroutine sim_outflowVelBlk2d(velOut,u,v,xcell,ycell,boundBox,dt,dx,dy,ix1,ix2,jy1,jy2)

   use Simulation_data
   use sim_outflowData

   implicit none
   real, intent(inout)                     :: velOut
   real, dimension(:,:,:), intent(in)      :: u,v
   real, dimension(:), intent(in)          :: xcell,ycell
   real, dimension(:,:), intent(in)        :: boundBox
   real, intent(in)                        :: dt,dx,dy
   integer, intent (in)                    :: ix1,ix2,jy1,jy2

   integer :: i,j,k
   real    :: xi,yi

   k = 1

   do j=jy1,jy2
    do i=ix1,ix2   
      xi = xcell(i)
      yi = ycell(j)

      if(yi .le. sim_yMax .and. yi .ge. sim_yMax-dy) velOut = max(velOut,v(i,j+1,k))

     end do
   end do
   
end subroutine sim_outflowVelBlk2d

subroutine sim_outflowVelBlk3d(velOut,u,v,w,xcell,ycell,zcell,boundBox,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)

   use Simulation_data
   use sim_outflowData

   implicit none
   real, intent(inout)                   :: velOut
   real, dimension(:,:,:), intent(in)    :: u,v,w
   real, dimension(:), intent(in)        :: xcell,ycell,zcell
   real, dimension(:,:), intent(in)      :: boundBox
   real, intent(in)                      :: dt,dx,dy,dz
   integer, intent (in)                  :: ix1,ix2,jy1,jy2,kz1,kz2

   integer :: i,j,k
   real    :: xi,yi,zi

   do k=kz1,kz2
     do j=jy1,jy2
       do i=ix1,ix2
        xi = xcell(i)
        yi = ycell(j)
        zi = zcell(k)

        if(yi .le. sim_yMax .and. yi .ge. sim_yMax-dy) velOut = max(velOut,v(i,j+1,k))

       end do
     end do
   end do

end subroutine sim_outflowVelBlk3d
