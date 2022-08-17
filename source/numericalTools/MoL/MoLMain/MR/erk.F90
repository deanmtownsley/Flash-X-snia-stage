module erk

   implicit none

contains

   subroutine euler_init()
      use mr_data, only: AF, bF, cF, mr_nstages_fast

      implicit none

      allocate (AF(1, 1))
      allocate (bF(1))
      allocate (cF(1))

      mr_nstages_fast = 1

      AF(1, 1) = 0d0

      bF(1) = 1d0

      cF(1) = 0d0
   end subroutine euler_init

   subroutine rk2_heun_init()
      use mr_data, only: AF, bF, cF, mr_nstages_fast

      implicit none

      allocate (AF(2, 2))
      allocate (bF(2))
      allocate (cF(2))

      mr_nstages_fast = 2

      AF = 0d0
      AF(2, 1) = 1d0

      bF(1) = 1d0/2d0
      bF(2) = 1d0/2d0

      cF(1) = 0d0
      cF(2) = 1d0
   end subroutine rk2_heun_init

   subroutine rk3_ssp_init()
      use mr_data, only: AF, bF, cF, mr_nstages_fast

      implicit none

      allocate (AF(3, 3))
      allocate (bF(3))
      allocate (cF(3))

      mr_nstages_fast = 3

      AF = 0d0
      AF(2, 1) = 1d0
      AF(3, 1) = 1d0/4d0
      AF(3, 2) = 1d0/4d0

      bF(1) = 1d0/6d0
      bF(2) = 1d0/6d0
      bF(3) = 2d0/3d0

      cF(1) = 0d0
      cF(2) = 1d0
      cF(3) = 1d0/2d0
   end subroutine rk3_ssp_init

   subroutine rk4_init()
      use mr_data, only: AF, bF, cF, mr_nstages_fast

      implicit none

      allocate (AF(4, 4))
      allocate (bF(4))
      allocate (cF(4))

      mr_nstages_fast = 4

      AF = 0d0
      AF(2, 1) = 1d0/2d0
      AF(3, 2) = 1d0/2d0
      AF(4, 3) = 1d0

      bF(1) = 1d0/6d0
      bF(2) = 1d0/3d0
      bF(3) = 1d0/3d0
      bF(4) = 1d0/6d0

      cF(1) = 0d0
      cF(2) = 1d0/2d0
      cF(3) = 1d0/2d0
      cF(4) = 1d0
   end subroutine rk4_init

   subroutine erk_init
      use mr_data, only: mr_fastMethod
      use ml_interface, only: ml_error

      implicit none

      select case (mr_fastMethod)

      case ("erk-euler")
         call euler_init

      case ("erk-rk2-heun")
         call rk2_heun_init

      case ("erk-rk3-ssp")
         call rk3_ssp_init

      case ("erk-rk4")
         call rk4_init

      case default
         call ml_error("Unkown ERK method")
      end select

   end subroutine erk_init

end module erk
