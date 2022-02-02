!!****if* source/Simulation/SimulationMain/incompFlow/PoolBoiling/localAPI/outflow/sim_outflowInterface
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
!! NAME
!!  
!!
!! SYNOPSIS
!!  sim_outflowInterface()
!!
!! DESCRIPTION
!!  This is an interface specific for outflow boundary conditions
!!
!!***

#include "constants.h"

Module sim_outflowInterface

  implicit none

  interface 
    subroutine sim_outflowInit()
    end subroutine sim_outflowInit
  end interface

  interface 
    subroutine sim_outflowFinalize()
    end subroutine sim_outflowFinalize
  end interface

  interface 
    subroutine sim_outflowSetBC(dt)
    real, intent(in) :: dt
    end subroutine sim_outflowSetBC
  end interface

  interface sim_outflowLSDampingBlk
    subroutine sim_outflowLSDampingBlk2d(pfrc,phi,xcell,ycell,boundBox,dt,dx,dy,ix1,ix2,jy1,jy2)
    real, dimension(:,:,:), intent(inout)  :: pfrc
    real, dimension(:,:,:), intent(in)     :: phi
    real, dimension(:), intent(in)         :: xcell,ycell
    real, dimension(:,:), intent(in)       :: boundBox
    real, intent(in)                       :: dt,dx,dy
    integer, intent(in)                    :: ix1,ix2,jy1,jy2
    end subroutine sim_outflowLSDampingBlk2d

    subroutine sim_outflowLSDampingBlk3d(pfrc,phi,xcell,ycell,zcell,boundBox,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)
    real, dimension(:,:,:), intent(inout)  :: pfrc
    real, dimension(:,:,:), intent(in)     :: phi
    real, dimension(:), intent(in)         :: xcell,ycell,zcell
    real, dimension(:,:), intent(in)       :: boundBox
    real, intent(in)                       :: dt,dx,dy,dz
    integer, intent (in)                   :: ix1,ix2,jy1,jy2,kz1,kz2
    end subroutine sim_outflowLSDampingBlk3d
  end interface

  interface sim_outflowVelBlk
    subroutine sim_outflowVelBlk2d(velOut,u,v,xcell,ycell,boundBox,dt,dx,dy,ix1,ix2,jy1,jy2)
    real, intent(inout)                    :: velOut
    real, dimension(:,:,:), intent(in)     :: u,v
    real, dimension(:), intent(in)         :: xcell,ycell
    real, dimension(:,:), intent(in)       :: boundBox
    real, intent(in)                       :: dt,dx,dy
    integer, intent(in)                    :: ix1,ix2,jy1,jy2
    end subroutine sim_outflowVelBlk2d

    subroutine sim_outflowVelBlk3d(velOut,u,v,w,xcell,ycell,zcell,boundBox,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)
    real, intent(inout)                    :: velOut
    real, dimension(:,:,:), intent(in)     :: u,v,w
    real, dimension(:), intent(in)         :: xcell,ycell,zcell
    real, dimension(:,:), intent(in)       :: boundBox
    real, intent(in)                       :: dt,dx,dy,dz
    integer, intent (in)                   :: ix1,ix2,jy1,jy2,kz1,kz2
    end subroutine sim_outflowVelBlk3d
  end interface

End module sim_outflowInterface
