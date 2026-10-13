!!****if* source/physics/ImBound/localAPI/ib_interface
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
!!***
module ib_interface

   implicit none

   !! This builds ANN tree for a given body data
   interface
      subroutine ib_annBuildTree(body)
         use ImBound_type, ONLY: ImBound_type_t
         type(ImBound_type_t), intent(INOUT)  :: body
      end subroutine ib_annBuildTree
   end interface

   interface
      subroutine ib_annDestroyTree(body)
         use ImBound_type, ONLY: ImBound_type_t
         type(ImBound_type_t), intent(INOUT)  :: body
      end subroutine ib_annDestroyTree
   end interface

   interface
      subroutine ib_annSearchTree(body, queryPt, annElems, annIdx)
         use ImBound_type, ONLY: ImBound_type_t
         type(ImBound_type_t), intent(IN)  :: body
         integer, intent(IN) :: annElems
         ! query point
         real, dimension(:), target, intent(IN) :: queryPt
         ! indices of nearest neighbors
         integer, dimension(:), target, intent(OUT):: annIdx
      end subroutine ib_annSearchTree
   end interface

   interface
      subroutine ib_readBody(body, bodyFile)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         type(ImBound_type_t), intent(inout) :: body
         character(len=*), intent(in) :: bodyFile
      end subroutine ib_readBody
   end interface

   interface ib_bruteForceMap
      subroutine ib_bruteForceMap2D(lmda, xcenter, ycenter, dx, dy, ix1, ix2, jy1, jy2, body)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         real, dimension(:, :, :), intent(inout) :: lmda
         real, dimension(:), intent(in) :: xcenter, ycenter
         type(ImBound_type_t), intent(in) :: body
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, intent(in) :: dx, dy
      end subroutine ib_bruteForceMap2D

      subroutine ib_bruteForceMap3D(lmda, xcenter, ycenter, zcenter, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, body)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         real, dimension(:, :, :), intent(inout) :: lmda
         real, dimension(:), intent(in) :: xcenter, ycenter, zcenter
         type(ImBound_type_t), intent(in) :: body
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: dx, dy, dz
      end subroutine ib_bruteForceMap3D
   end interface ib_bruteForceMap

   interface ib_annMap
      subroutine ib_annMap2D(lmda, xcenter, ycenter, dx, dy, ix1, ix2, jy1, jy2, body)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         real, dimension(:, :, :), intent(inout) :: lmda
         real, dimension(:), intent(in) :: xcenter, ycenter
         type(ImBound_type_t), intent(in) :: body
         integer, intent(in) :: ix1, ix2, jy1, jy2
         real, intent(in) :: dx, dy
      end subroutine ib_annMap2D

      subroutine ib_annMap3D(lmda, xcenter, ycenter, zcenter, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, body)
         use ImBound_type, ONLY: ImBound_type_t
         implicit none
         real, dimension(:, :, :), intent(inout) :: lmda
         real, dimension(:), intent(in) :: xcenter, ycenter, zcenter
         type(ImBound_type_t), intent(in) :: body
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
         real, intent(in) :: dx, dy, dz
      end subroutine ib_annMap3D
   end interface

   interface
      subroutine ib_velGfm2d_fixed(lmda, velx, vely, px, py, dt, coeff, buffer, &
                                   dx, dy, ix1, ix2, jy1, jy2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: velx, vely
         real, dimension(:, :, :), intent(in) :: lmda
         real, dimension(:, :, :), intent(in) :: px, py
         real, intent(in) :: dt, dx, dy, coeff, buffer(3)
         integer, intent(in) :: ix1, ix2, jy1, jy2
      end subroutine ib_velGfm2d_fixed
   end interface

   interface
      subroutine ib_velGfm3d_fixed(lmda, velx, vely, velz, px, py, pz, &
                                   dt, coeff, buffer, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2)
         implicit none
         real, dimension(:, :, :), intent(inout) :: velx, vely, velz
         real, dimension(:, :, :), intent(in) :: lmda
         real, dimension(:, :, :), intent(in) :: px, py, pz
         real, intent(in) :: dt, dx, dy, dz, coeff, buffer(3)
         integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
      end subroutine ib_velGfm3d_fixed
   end interface

end module ib_interface
