!! source/physics/ImBound/localAPI/ib_annMap
!!
!! NAME
!!
!! ib_annMap(blockCount,blockList,dt)
!!
!! SYNOPSIS
!!
!!
!! VARIABLES
!!
!!
!! DESCRIPTION
!!
!! Subroutine to find the distance function lambda for
!! the immersed boundary (IB).
!!
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
