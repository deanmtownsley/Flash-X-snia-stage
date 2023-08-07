!! source/physics/ImBound/ImBoundMain/ib_annMap
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
subroutine ib_annMap(lmda, xcenter, ycenter, dx, dy, ix1, ix2, jy1, jy2, body)

#include "Simulation.h"
#include "constants.h"

   ! Modules Used
   use ImBound_data, ONLY: ib_annQueries
   use ImBound_type, ONLY: ImBound_type_t
   use ib_interface, ONLY: ib_annSearchTree
   use Timers_interface, ONLY: Timers_start, Timers_stop
   implicit none

   ! Arguments
   real, dimension(:, :, :), intent(inout) :: lmda
   real, dimension(:), intent(in) :: xcenter, ycenter
   type(ImBound_type_t), intent(in) :: body
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: dx, dy

   ! Internal variables
   integer :: i, j, k, panelIndex
   real :: xcell, ycell, zcell, mvd

   ! For the algorithm
   real, allocatable, dimension(:) :: PA, PB, Pcell, P0, v1
   real, allocatable, dimension(:) :: dist
   real :: u

   ! the grid/search data
   integer :: annElems, annIndex
   ! indices of nearest neighbors (NN)
   integer, dimension(:), allocatable :: annIdx
   real :: eps = 1e-13, dotNorm, magNorm, thetaNorm

   ! allocating data
   allocate (PA(body%dims), PB(body%dims), Pcell(body%dims), P0(body%dims), v1(body%dims))

   ! define ANN parameters
   annElems = ib_annQueries
   allocate (dist(annElems))
   allocate (annIdx(annElems))

   k = 1
   do j = jy1, jy2
      do i = ix1, ix2
         dist = 1e+200
         mvd = 1e+200
         PA = 0.0
         PB = 0.0
         P0 = 0.0
         Pcell = 0.0

         v1 = 0.0

         ! x and y coordinates for the current grid cell
         xcell = xcenter(i)
         ycell = ycenter(j)
         zcell = 0.0

         ! Grid cell point
         Pcell = (/xcell, ycell/)

         !!! find the  nearest neighbors
         !! to compute ls value
         call Timers_start("ib_annSearchTree")
         call ib_annSearchTree(body, Pcell, annElems, annIdx)
         call Timers_stop("ib_annSearchTree")

         do annIndex = 1, annElems
            panelIndex = annIdx(annIndex) + 1 ! need + 1 to convert c++ index to fortran

            PA = body%elems(panelIndex)%pA
            PB = body%elems(panelIndex)%pB

            ! Drop a normal from Pcell to the line made by connecting PA PB (not the
            ! line segment)
            u = ((Pcell(1) - PA(1))*(PB(1) - PA(1)) + (Pcell(2) - PA(2))*(PB(2) - PA(2)))/ &
                (((PB(1) - PA(1))**2) + ((PB(2) - PA(2))**2))

            ! Re-assign u if the normal hits the line segment to the left of PA or
            ! the right of PB
            if (u .lt. 0) then
               u = 0.0
            else if (u .gt. 1) then
               u = 1.0
            end if

            ! Find the point on the line segment with the shortest distance to Pcell
            ! (If the normal hits the line outside the line segment it is
            !  reassigned to hit the closer endpoint.)
            P0(1) = PA(1) + (PB(1) - PA(1))*u
            P0(2) = PA(2) + (PB(2) - PA(2))*u

            ! Determine the quadrent and angle for the "normal"
            ! (If to the left or right of the line segment the vector with the
            !  shortest distance to the line segment will not be perpendicular)

            if (abs(P0(1) - PA(1)) .lt. 1e-13 .and. abs(P0(2) - PA(2)) .lt. 1e-13) then
               v1 = (/(Pcell(1) - P0(1)), (Pcell(2) - P0(2))/)
            else
               v1 = (/(Pcell(1) - P0(1)), (Pcell(2) - P0(2))/)
            end if

            dotNorm = v1(1)*body%elems(panelIndex)%normal(1) + v1(2)*body%elems(panelIndex)%normal(2)
            magNorm = sqrt(body%elems(panelIndex)%normal(1)**2 + body%elems(panelIndex)%normal(2)**2)* &
                      sqrt(v1(1)**2 + v1(2)**2)

            thetaNorm = acos(dotNorm/(magNorm + eps))

            if ((thetaNorm <= (acos(-1.)/2 + eps)) .or. (thetaNorm > (3*acos(-1.)/2 - eps))) then
               dist(annIndex) = -sqrt(v1(1)**2 + v1(2)**2)
            else
               dist(annIndex) = sqrt(v1(1)**2 + v1(2)**2)
            end if

            if (abs(mvd) > abs(dist(annIndex))) then
               mvd = dist(annIndex)
            end if

         end do

         ! For first body explicitly satisfy level set, and then compare with
         ! existing level set for successive bodies
         lmda(i, j, k) = mvd

      end do
   end do

   deallocate (dist)
   deallocate (PA, PB, Pcell, P0, v1)
   deallocate (annIdx)

end subroutine ib_annMap
