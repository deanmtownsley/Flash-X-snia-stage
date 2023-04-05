!!****if* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/rt_tm_interface
!!
!! NAME
!!   rt_tm_interface
!!
!! SYNOPSIS
!!   use rt_tm_interface : ONLY
!!
!!  DESCRIPTION 
!!    Interface for internal RadTrans TwoMoment subroutines
!!
!!***
module rt_tm_interface

  implicit none

#include "constants.h"

  interface
     subroutine rt_tm_reconstruction(Uin,nX,lo,hi,u_lo,u_hi)
        real, pointer, dimension(:,:,:,:) :: Uin
        integer, dimension(3), intent(in) :: nX
        integer, dimension(MDIM), intent(in) :: lo,hi,u_lo,u_hi
     end subroutine rt_tm_reconstruction
  end interface

  interface
     subroutine rt_tm_projection(Uin,Sout,nX,lo,hi,u_lo,u_hi)
        real, pointer, dimension(:,:,:,:) :: Uin,Sout
        integer, dimension(3), intent(in) :: nX
        integer, dimension(MDIM), intent(in) :: lo,hi,u_lo,u_hi
     end subroutine rt_tm_projection
  end interface

end module rt_tm_interface
