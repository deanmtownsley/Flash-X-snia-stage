!!****if* source/physics/sourceTerms/Flame/FlameMain/RDSplit5point/fl_laplacian
!!
!! NAME
!!
!!  fl_laplacian
!!
!! SYNOPSIS
!!
!!  call fl_laplacian(real, dimension(:,:,:)(out) :: lapl,
!!                    real, dimension(:,:,:)(in) :: flam,
!!                    integer(in) :: h,
!!                    type(Grid_tile_t) :: tileDesc)
!!
!! DESCRIPTION
!!
!! calculate the laplacian of the flame progress variable
!!
!! Dean Townsley 2008
!!
!!  The block id (bid) is passed in so that we can retrieve coordinate info
!!  The laplacian is only computed for the interior cells, although the
!!  indices of lapl also run over the guard cells to simplify indexing
!!
!! ARGUMENTS
!!
!!   lapl : laplacian 
!!
!!   flam : flame 
!!
!!   h : step size
!!
!!   bid : block id
!!
!!
!!
!!***


#include "Simulation.h"
#include "constants.h"
subroutine fl_laplacian(lapl, flam, h, tileDesc)

  use Grid_tile, only: Grid_tile_t
  use Grid_interface, only : Grid_getGeometry, Grid_getBlkIndexLimits, &
                             Grid_getDeltas, Grid_getCellCoords
  use Driver_interface, only : Driver_abort
  implicit none
  real, dimension(:,:,:), intent(out) :: lapl
  real, dimension(:,:,:), intent(in) :: flam
  integer, intent(in) :: h
  type(Grid_tile_t) :: tileDesc

  integer :: geom
  real, dimension(MDIM) :: celldeltas, deltas

  real :: inv_12_dx2, inv_12_dy2, inv_12_dz2
  real :: inv_12_dr2, inv_12_dr, inv_12_dtheta, inv_12_dtheta2, inv_12_dphi2
  real :: d2, d1
  real, dimension(:), allocatable :: inv_r, inv_r2, two_over_r, ctan, inv_sin2
  integer :: i,j,k, istat
  integer :: h2

  h2 = 2*h

  if (h2 > NGUARD) call Driver_abort("Step size in fl_laplacian is too large for the number of guard cells")

  call Grid_getGeometry(geom)
  tileDesc%deltas(celldeltas)
  deltas(IAXIS) = h * celldeltas(IAXIS)
  if (NDIM >= 2) deltas(JAXIS) = h * celldeltas(JAXIS)
  if (NDIM == 3) deltas(KAXIS) = h * celldeltas(KAXIS)

  select case (geom)
  case (CARTESIAN)

     ! seems like any self-respecting optimizing compiler would be able to
     ! pull these out of the loop but we won't trust that
     inv_12_dx2 = 1.0/12.0/deltas(IAXIS)**2
     if (NDIM >= 2) inv_12_dy2 = 1.0/12.0/deltas(JAXIS)**2
     if (NDIM == 3) inv_12_dz2 = 1.0/12.0/deltas(KAXIS)**2

     do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
           
              lapl(i,j,k) = ( -flam(i-h2,j,k) + 16*flam(i-h,j,k) -30*flam(i,j,k) &
                                             + 16*flam(i+h,j,k) - flam(i+h2,j,k) ) * inv_12_dx2
              if (NDIM >= 2) then
                 lapl(i,j,k) = lapl(i,j,k) + &
                            ( -flam(i,j-h2,k) + 16*flam(i,j-h,k) -30*flam(i,j,k) &
                                             + 16*flam(i,j+h,k) - flam(i,j+h2,k) ) * inv_12_dy2
              endif
              if (NDIM == 3) then
                 lapl(i,j,k) = lapl(i,j,k) + &
                            ( -flam(i,j,k-h2) + 16*flam(i,j,k-h) -30*flam(i,j,k) &
                                             + 16*flam(i,j,k+h) - flam(i,j,k+h2) ) * inv_12_dz2
              endif

           enddo
        enddo
     enddo

  case (CYLINDRICAL)

     inv_12_dr  = 1.0/12.0/deltas(IAXIS)
     inv_12_dr2 = inv_12_dr/deltas(IAXIS)

     allocate( inv_r( tileDesc%limits(LOW,IAXIS) : tileDesc%limits(HIGH,IAXIS) ), STAT=istat )
     if (istat/=0) call Driver_abort("Unable to allocate inv_r in fl_laplacian")
     do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
        inv_r(i) = 1.0/celldeltas(IAXIS)/(i-0.5)
     enddo

     if (NDIM >= 2) inv_12_dz2 = 1.0/12.0/deltas(JAXIS)**2
     if (NDIM == 3) then
        inv_12_dtheta2 = 1.0/12.0/deltas(KAXIS)**2
        allocate( inv_r2( tileDesc%limits(LOW,IAXIS) : tileDesc%limits(HIGH,IAXIS) ), STAT=istat )
        if (istat/=0) call Driver_abort("Unable to allocate inv_r2 in fl_laplacian")
        inv_r2(:) = inv_r(:)**2
     endif


     do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)

              d2 = ( -flam(i-h2,j,k) + 16*flam(i-h,j,k) -30*flam(i,j,k) &
                                    + 16*flam(i+h,j,k) - flam(i+h2,j,k) ) * inv_12_dr2
              d1 = ( flam(i-h2,j,k) -8*flam(i-h,j,k)  +8*flam(i+h,j,k) -flam(i+h2,j,k) )*inv_12_dr

              lapl(i,j,k) = d2 + inv_r(i)*d1
                            
              if (NDIM >= 2) then
                 d2 = ( -flam(i,j-h2,k) + 16*flam(i,j-h,k) -30*flam(i,j,k) &
                                       + 16*flam(i,j+h,k) - flam(i,j+h2,k) ) * inv_12_dz2
                 lapl(i,j,k) = lapl(i,j,k) +  d2
              endif
              if (NDIM == 3) then
                 d2 = ( -flam(i,j,k-h2) + 16*flam(i,j,k-h) -30*flam(i,j,k) &
                                       + 16*flam(i,j,k+h) - flam(i,j,k+h2) ) * inv_12_dtheta2
                 lapl(i,j,k) = lapl(i,j,k) +  inv_r2(i)*d2
              endif

           enddo
        enddo
     enddo

     if (NDIM==3) deallocate(inv_r2)
     deallocate(inv_r)

  case (SPHERICAL)

     inv_12_dr  = 1.0/12.0/deltas(IAXIS)
     inv_12_dr2 = inv_12_dr/deltas(IAXIS)

     allocate( inv_r2( tileDesc%limits(LOW:IAXIS) : tileDesc%limits(HIGH,IAXIS) ), STAT=istat )
     if (istat/=0) call Driver_abort("Unable to allocate inv_r2 in fl_laplacian")
     allocate( two_over_r( tileDesc%limits(LOW:IAXIS) : tileDesc%limits(HIGH,IAXIS) ), STAT=istat )
     if (istat/=0) call Driver_abort("Unable to allocate two_over_r in fl_laplacian")
     call Grid_getCellCoords(IAXIS, bid, CENTER, .true., inv_r2, blkLimitsGC(HIGH,IAXIS))
     do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
        inv_r2(i) = 1.0/celldeltas(IAXIS)/(i-0.5)
        two_over_r(i) = 2.0*inv_r2(i)
        inv_r2(i) = inv_r2(i)**2
     enddo

     if (NDIM >= 2) then
        inv_12_dtheta = 1.0/12.0/deltas(JAXIS)
        inv_12_dtheta2 = inv_12_dtheta/deltas(JAXIS)
        allocate( ctan( tileDesc%limits(LOW:JAXIS) : tileDesc%limits(HIGH,JAXIS) ), STAT=istat )
        if (istat/=0) call Driver_abort("Unable to allocate ctan in fl_laplacian")
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           ctan(j) = 1.0/tan( celldeltas(JAXIS)*(j-0.5))
        enddo
     endif 

     if (NDIM == 3) then
        inv_12_dphi2 = 1.0/12.0/deltas(KAXIS)**2
        allocate( inv_sin2( tileDesc%limits(LOW:JAXIS) : tileDesc%limits(HIGH,JAXIS) ), STAT=istat )
        if (istat/=0) call Driver_abort("Unable to allocate inv_sin2 in fl_laplacian")
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           inv_sin2(j) = 1.0+ctan(j)**2
        enddo
     endif


     do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
                d2 = ( -flam(i-h2,j,k) + 16*flam(i-h,j,k) -30*flam(i,j,k) &
                                      + 16*flam(i+h,j,k) - flam(i+h2,j,k) ) * inv_12_dr2
                d1 = ( flam(i-h2,j,k)-8*flam(i-h,j,k)+8*flam(i+h,j,k) -flam(i+h2,j,k) )*inv_12_dr
                lapl(i,j,k) = d2 + two_over_r(i)*d1

                if (NDIM>=2) then
                   d2 = ( -flam(i,j-h2,k) + 16*flam(i,j-h,k) -30*flam(i,j,k) &
                                         + 16*flam(i,j+h,k) - flam(i,j+h2,k) ) * inv_12_dtheta2
                   d1 = ( flam(i,j-h2,k)-8*flam(i,j-h,k)+8*flam(i,j+h,k)-flam(i,j+h2,k) )*inv_12_dtheta
                   lapl(i,j,k) = lapl(i,j,k) + inv_r2(i)*( d2 + ctan(j)*d1 )
                endif
                if (NDIM==3) then
                   d2 = ( -flam(i,j,k-h2) + 16*flam(i,j,k-h) -30*flam(i,j,k) &
                                         + 16*flam(i,j,k+h) - flam(i,j,k+h2) ) * inv_12_dphi2
                   lapl(i,j,k) = lapl(i,j,k) + inv_r2(i)*inv_sin2(j)*d2
                endif
             enddo
          enddo
       enddo

       if (NDIM==3) deallocate(inv_sin2)
       if (NDIM>=2) deallocate(ctan)
       deallocate(inv_r2)
       deallocate(two_over_r)

 
  end select

  return
end subroutine
