module vector

   implicit none

   private
   public :: vec_magnitude2D, vec_magnitude3D, vec_cross_product3D
   public :: vec_rotate2D, vec_offset2D, vec_offset3D

contains

   function vec_magnitude2D(vec) result(mag)
      implicit none
      real, dimension(2), intent(in)  :: vec
      real :: mag
      mag = sqrt(vec(1)**2 + vec(2)**2)
   end function vec_magnitude2D

   function vec_magnitude3D(vec) result(mag)
      implicit none
      real, dimension(3), intent(in)  :: vec
      real :: mag
      mag = sqrt(vec(1)**2 + vec(2)**2 + vec(3)**2)
   end function vec_magnitude3D

   function vec_cross_product3D(a, b) result(cross)
      implicit none
      real, dimension(3), intent(in) :: a, b
      real, dimension(3) :: cross
      cross(1) = a(2)*b(3) - a(3)*b(2)
      cross(2) = a(3)*b(1) - a(1)*b(3)
      cross(3) = a(1)*b(2) - a(2)*b(1)
   end function vec_cross_product3D

   subroutine vec_rotate2D(vec, rotate)
      implicit none
      real, dimension(2), intent(inout) :: vec
      real, dimension(2, 2), intent(in) :: rotate
      vec = matmul(vec, rotate)
   end subroutine vec_rotate2D

   subroutine vec_offset2D(vec, offset)
      implicit none
      real, dimension(2), intent(inout) :: vec
      real, dimension(2), intent(in) :: offset
      vec = vec + offset
   end subroutine vec_offset2D

   subroutine vec_offset3D(vec, offset)
      implicit none
      real, dimension(3), intent(inout) :: vec
      real, dimension(3), intent(in) :: offset
      vec = vec + offset
   end subroutine vec_offset3D

end module
