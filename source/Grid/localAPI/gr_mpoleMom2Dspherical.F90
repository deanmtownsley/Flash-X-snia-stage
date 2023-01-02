!!****if* source/Grid/GridSolvers/Multipole_new/gr_mpoleMom2Dspherical
!!
!! NAME
!!
!!  gr_mpoleMom2Dspherical
!!
!! SYNOPSIS
!!
!!  gr_mpoleMom2Dspherical (integer (in) :: idensvar)
!!
!! DESCRIPTION
!!
!!  Prepares for evaluation of the moments in 2D spherical geometry. In this
!!  routine, all the necessary arrays are prepared to enable evaluation of
!!  the moments in radial bin order. Each of the moments are grouped together
!!  according to their radial bins. This will ensure optimum unit stride values
!!  when accessing the big moment arrays and makes threading trivial.
!!
!! ARGUMENTS
!!
!!  idensvar : the index of the density variable
!!
!!***

!!REORDER(4): solnData

subroutine gr_mpoleMom2Dspherical (idensvar)

  implicit none
  
  integer, intent (in) :: idensvar

  return
end subroutine gr_mpoleMom2Dspherical
