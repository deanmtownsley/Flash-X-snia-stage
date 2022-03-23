!!****h* source/numericalTools/Stencils/Stencils_interface
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
!!  Stencils_interface
!!
!! SYNOPSIS
!!
!!   use Stencils_interface
!!
!! DESCRIPTION
!!
!!  This is the module needed for defining the Stencils unit interface.
!!
!!***

Module Stencils_interface

  interface
    subroutine Stencils_finalize()
    end subroutine Stencils_finalize
  end interface

  interface
    subroutine Stencils_init()
    end subroutine Stencils_init
  end interface

  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  interface Stencils_integrateEuler
     subroutine Stencils_integrateEulerScalar(phi,rhs,dt,ix1,ix2,jy1,jy2,kz1,kz2,iSource)
     implicit none
     real, dimension(:,:,:), intent(inout):: phi
     real, dimension(:,:,:), intent(in) :: rhs
     real, intent(in) :: dt
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
     real, intent(in) :: iSource
     end subroutine Stencils_integrateEulerScalar

     subroutine Stencils_integrateEulerArray(phi,rhs,dt,ix1,ix2,jy1,jy2,kz1,kz2,iSource)
     implicit none
     real, dimension(:,:,:), intent(inout):: phi
     real, dimension(:,:,:), intent(in) :: rhs
     real, intent(in) :: dt
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
     real, dimension(:,:,:), intent(in) :: iSource
     end subroutine Stencils_integrateEulerArray
  end interface

  interface Stencils_integrateAB2
     subroutine Stencils_integrateAB2Scalar(phi,rhsNew,rhsOld,dt,ix1,ix2,jy1,jy2,kz1,kz2,iSource)
     implicit none
     real, dimension(:,:,:), intent(inout):: phi
     real, dimension(:,:,:), intent(in) :: rhsNew, rhsOld
     real, intent(in) :: dt
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
     real, intent(in) :: iSource
     end subroutine Stencils_integrateAB2Scalar

     subroutine Stencils_integrateAB2Array(phi,rhsNew,rhsOld,dt,ix1,ix2,jy1,jy2,kz1,kz2,iSource)
     implicit none
     real, dimension(:,:,:), intent(inout):: phi
     real, dimension(:,:,:), intent(in) :: rhsNew, rhsOld
     real, intent(in) :: dt
     integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
     real, dimension(:,:,:), intent(in) :: iSource
     end subroutine Stencils_integrateAB2Array
  end interface

  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  interface
    subroutine Stencils_cnt_advectUpwind2d(rhs,phi,u,v,dx,dy,ix1,ix2,jy1,jy2)
    implicit none
    real, dimension(:,:,:), intent(inout):: rhs
    real, dimension(:,:,:), intent(in) :: phi,u,v
    real, intent(in) :: dx,dy
    integer, intent(in) :: ix1,ix2,jy1,jy2
    end subroutine Stencils_cnt_advectUpwind2d
  end interface

  interface
    subroutine Stencils_cnt_advectUpwind3d(rhs,phi,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)
    implicit none
    real, dimension(:,:,:), intent(inout):: rhs
    real, dimension(:,:,:), intent(in) :: phi,u,v,w
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    end subroutine Stencils_cnt_advectUpwind3d
  end interface

  interface
    subroutine Stencils_advectWeno2d(rhs,phi,u,v,dx,dy,ix1,ix2,jy1,jy2,center,facex,facey)
    implicit none
    real, dimension(:,:,:), intent(inout):: rhs
    real, dimension(:,:,:), intent(in) :: phi,u,v
    real, intent(in) :: dx,dy
    integer, intent(in) :: ix1,ix2,jy1,jy2
    integer, intent(in) :: center,facex,facey
    end subroutine Stencils_advectWeno2d
  end interface
 
  interface
    subroutine Stencils_advectWeno3d(rhs,phi,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,&
                                      center,facex,facey,facez)
    implicit none
    real, dimension(:,:,:), intent(inout):: rhs
    real, dimension(:,:,:), intent(in) :: phi,u,v,w
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    integer, intent(in) :: center,facex,facey,facez
    end subroutine Stencils_advectWeno3d
  end interface

  interface
    subroutine Stencils_advectCentral2d(rhs,phi,u,v,dx,dy,ix1,ix2,jy1,jy2,center,facex,facey)
    implicit none
    real, dimension(:,:,:), intent(inout):: rhs
    real, dimension(:,:,:), intent(in) :: phi,u,v
    real, intent(in) :: dx,dy
    integer, intent(in) :: ix1,ix2,jy1,jy2
    integer, intent(in) :: center,facex,facey
    end subroutine Stencils_advectCentral2d
  end interface

  interface
    subroutine Stencils_advectCentral3d(rhs,phi,u,v,w,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2,&
                                         center,facex,facey,facez)
    implicit none
    real, dimension(:,:,:), intent(inout):: rhs
    real, dimension(:,:,:), intent(in) :: phi,u,v,w
    real, intent(in) :: dx,dy,dz
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    integer, intent(in) :: center,facex,facey,facez
    end subroutine Stencils_advectCentral3d
  end interface

  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  interface Stencils_diffusion2d
    subroutine Stencils_diffusion2dConst(rhs, phi, dx, dy, Coeff, ix1, ix2, jy1, jy2)
    implicit none
    real, dimension(:,:,:), intent(inout) :: rhs
    real, dimension(:,:,:), intent(in)  :: phi
    real, intent(in) :: Coeff
    real, intent(in) :: dx, dy
    integer, intent(in) :: ix1, ix2, jy1, jy2
    end subroutine Stencils_diffusion2dConst

    subroutine Stencils_diffusion2dVar(rhs, phi, dx, dy, Coeff, ix1, ix2, jy1, jy2)
    implicit none
    real, dimension(:,:,:), intent(inout) :: rhs
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(in)  :: Coeff
    real, intent(in) :: dx, dy
    integer, intent(in) :: ix1, ix2, jy1, jy2
    end subroutine Stencils_diffusion2dVar
  end interface

  interface Stencils_diffusion3d
    subroutine Stencils_diffusion3dConst(rhs, phi, dx, dy, dz, Coeff, ix1, ix2, jy1, jy2, kz1, kz2)
    implicit none
    real, dimension(:,:,:), intent(inout) :: rhs
    real, dimension(:,:,:), intent(in)  :: phi
    real, intent(in) :: dx, dy, dz
    real, intent(in) :: Coeff
    integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
    end subroutine Stencils_diffusion3dConst

    subroutine Stencils_diffusion3dVar(rhs, phi, dx, dy, dz, Coeff, ix1, ix2, jy1, jy2, kz1, kz2)
    implicit none
    real, dimension(:,:,:), intent(inout) :: rhs
    real, dimension(:,:,:), intent(in)  :: phi
    real, intent(in) :: dx, dy, dz
    real, dimension(:,:,:), intent(in) :: Coeff
    integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
    end subroutine Stencils_diffusion3dVar
  end interface

  !----------------------------------------------------------------------------------------------------------
  !----------------------------------------------------------------------------------------------------------

  interface
    subroutine Stencils_lsRedistance2d(phi,phi_orig,dt,dx,dy,ix1,ix2,jy1,jy2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, dimension(:,:,:), intent(inout):: phi
    real, intent(in) :: dx,dy,dt
    real, dimension(:,:,:), intent(in):: phi_orig
    end subroutine Stencils_lsRedistance2d
  end interface

  interface
    subroutine Stencils_lsRedistance3d(phi,phi_orig,dt,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, dimension(:,:,:), intent(inout):: phi
    real, intent(in) :: dx,dy,dz,dt
    real, dimension(:,:,:), intent(in):: phi_orig
    end subroutine Stencils_lsRedistance3d
  end interface

  interface Stencils_lsFaceProps2d
    subroutine Stencils_lsFaceProps2dFixedSmear(phi,iPropX,iPropY,lsScalarProp,ix1,ix2,jy1,jy2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY
    end subroutine Stencils_lsFaceProps2dFixedSmear

    subroutine Stencils_lsFaceProps2dUserSmear(phi,iPropX,iPropY,lsScalarProp,ix1,ix2,jy1,jy2,iSmear)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY
    real, intent(in) :: iSmear
    end subroutine Stencils_lsFaceProps2dUserSmear
  end interface

  interface Stencils_lsFaceProps3d
    subroutine Stencils_lsFaceProps3dFixedSmear(phi,iPropX,iPropY,iPropZ,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY,iPropZ
    end subroutine Stencils_lsFaceProps3dFixedSmear

    subroutine Stencils_lsFaceProps3dUserSmear(phi,iPropX,iPropY,iPropZ,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2,iSmear)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropX,iPropY,iPropZ
    real, intent(in) :: iSmear
    end subroutine Stencils_lsFaceProps3dUserSmear
  end interface

  interface Stencils_lsCenterProps
    subroutine Stencils_lsCenterPropsSharp(phi,iPropC,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropC
    end subroutine Stencils_lsCenterPropsSharp

    subroutine Stencils_lsCenterPropsSmeared(phi,iPropC,lsScalarProp,ix1,ix2,jy1,jy2,kz1,kz2,iSmear)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, intent(in) :: lsScalarProp
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iPropC
    real, intent(in) :: iSmear
    end subroutine Stencils_lsCenterPropsSmeared
  end interface

  interface
    subroutine Stencils_lsNormals2d(phi,iNrmx,iNrmy,dx,dy,ix1,ix2,jy1,jy2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iNrmx, iNrmy
    real, intent(in) :: dx,dy
    end subroutine Stencils_lsNormals2d
  end interface

  interface
    subroutine Stencils_lsNormals3d(phi,iNrmx,iNrmy,iNrmz,dx,dy,dz,ix1,ix2,jy1,jy2,kz1,kz2)
    implicit none
    integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
    real, dimension(:,:,:), intent(in)  :: phi
    real, dimension(:,:,:), intent(inout) :: iNrmx, iNrmy, iNrmz
    real, intent(in) :: dx,dy,dz
    end subroutine Stencils_lsNormals3d
  end interface

end Module Stencils_interface
