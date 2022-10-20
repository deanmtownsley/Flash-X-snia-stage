!! source/physics/ImBound/ImBoundMain/ib_mapToGrid
!!
!! NAME
!!
!! ib_mapToGrid(blockCount,blockList,dt)
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
subroutine ib_mapToGrid(lmda, xcenter, ycenter, dx, dy, ix1, ix2, jy1, jy2, body)

#include "Simulation.h"
#include "constants.h"

   ! Modules Used
   use ImBound_type, ONLY: ImBound_type_t
   implicit none

   ! Arguments
   real, dimension(:, :, :), intent(inout) :: lmda
   real, dimension(:), intent(in) :: xcenter, ycenter
   type(ImBound_type_t), intent(in) :: body
   integer, intent(in) :: ix1, ix2, jy1, jy2
   real, intent(in) :: dx, dy

   ! Internal Variables
   integer :: numPart, e, ptelem, nel, p
   integer, allocatable, dimension(:) :: max_ptelem

   integer :: i, j, k, p_i, ap_i
   real :: mva, mvd
   real :: xcell, ycell, zcell

   ! For the algorithm
   real, allocatable, dimension(:) :: PA, PB, Pcell, P0, v1, v2
   real, allocatable, dimension(:) :: dist, distls
   real :: u
   integer :: nelm = 2 ! Dimension for the points, 2 for (x,y) in 2-D
   integer :: countit
   real    :: miny, maxy, mratio, nratio, xit

   ! the grid/search data
   integer :: nn_i, nn, nn_irc, nn_rc
   ! query points
   real, dimension(:), allocatable:: queryPt, queryPt_rc
   ! distance of queryPt from neighbors
   real, dimension(:), allocatable :: dists, dists_rc
   ! indices of nearest neighbors (NN)
   integer, dimension(:), allocatable :: nnIdx, nnIdx_rc
   real :: eps
   ! allocating data
   allocate (dist(body%numElems))
   allocate (PA(nelm), PB(nelm), Pcell(nelm), P0(nelm), v1(nelm), v2(nelm))
   !! define ANN parameters
   eps = 0.0
   nn = 3 ! #NN for dist value
   nn_rc = 8 ! #NN for ray casting
   allocate (distls(nn))
   !! allocate query point
   allocate (queryPt(body%dim))
   allocate (queryPt_rc(1))
   k = 1
   do j = jy1, jy2
      do i = ix1, ix2
         dist = 1e+200
         distls = 1e+200
         countit = 0   ! Counter to check no. of intersections with the body
         PA = 0.0
         PB = 0.0
         P0 = 0.0
         Pcell = 0.0

         v1 = 0.0
         v2 = 0.0

         ! x and y coordinates for the current grid cell
         xcell = xcenter(i)
         ycell = ycenter(j)
         zcell = 0.0
         !! ======================================= ANN approach ==============================================
         if (.true.) then
            queryPt = (/xcell, ycell/)
            queryPt_rc = ycell
            !!! find the  nearest neighbors
            !! to compute ls value
            call body%searchTree(queryPt, nn, nnIdx, dists, eps)
            !! for ray casting
            call body%searchTreeRC(queryPt_rc, nn_rc, nnIdx_rc, dists_rc, eps)
            do nn_i = 1, nn
               p_i = nnIdx(nn_i) + 1 ! need + 1 to convert c++ index to fortran
               PA = (/body%elems(p_i)%xA, body%elems(p_i)%yA/)
               PB = (/body%elems(p_i)%xB, body%elems(p_i)%yB/)
               ! Grid cell point
               Pcell = (/xcell, ycell/)

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
                  v2 = (/(P0(1) - PB(1)), (P0(2) - PB(2))/)
               else
                  v1 = (/(Pcell(1) - P0(1)), (Pcell(2) - P0(2))/)
                  v2 = (/(PA(1) - P0(1)), (PA(2) - P0(2))/)
               end if

               distls(nn_i) = sqrt(v1(1)**2 + v1(2)**2)

               !!! ray casting
               do nn_irc = 1, nn_rc
                  ! End points for the line segment of the IB
                  ! PA is on the left and PB is on the right
                  ap_i = nnIdx_rc(nn_irc) + 1 !! ap_i is short for panel_index
                  PA = (/body%elems(ap_i)%xA, body%elems(ap_i)%yA/)
                  PB = (/body%elems(ap_i)%xB, body%elems(ap_i)%yB/)
                  ! Grid cell point
                  Pcell = (/xcell, ycell/)
                  ! Find if the horizontal ray on right-side intersects with body
                  miny = min(PA(2), PB(2))
                  maxy = max(PA(2), PB(2))

                  if (ycell .gt. miny .and. ycell .lt. maxy) then

                     ! Method #1 use ratios to divide the current panel using
                     ! y intersection and find x

                     mratio = PA(2) - ycell
                     nratio = ycell - PB(2)
                     xit = (mratio*PB(1) + nratio*PA(1))/(mratio + nratio)

                     ! Method #2 use the equation of line instead

                     !mratio = (PB(2)-PA(2))/(PB(1)-PA(1))
                     !xit = PA(1) + (ycell - PA(2))/mratio

                     ! Check to make sure that the intersection is on the right

                     if (xit .ge. xcell) countit = countit + 1

                  end if
               end do
            end do
            ! Construct level set - if intersections are positive then the point
            ! lies outside (-), if odd then the point lies inside (+)
            mvd = sign(minval(distls(:)), 2*mod(countit, 2) - 1.)
            ! For first body explicitly satisfy level set, and then compare with
            ! existing level set for successive bodies
            lmda(i, j, k) = mvd
         end if
         !! ======================================= ANN approach ends ===============================================

         !! ======================================= Classical approach ==============================================
         if (.false.) then
            do p_i = 1, body%numElems ! p_i is short for panel_index
               ! End points for the line segment of the IB
               ! PA is on the left and PB is on the right
               PA = (/body%elems(p_i)%xA, body%elems(p_i)%yA/)
               PB = (/body%elems(p_i)%xB, body%elems(p_i)%yB/)
               ! Grid cell point
               Pcell = (/xcell, ycell/)
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
                  v2 = (/(P0(1) - PB(1)), (P0(2) - PB(2))/)
               else
                  v1 = (/(Pcell(1) - P0(1)), (Pcell(2) - P0(2))/)
                  v2 = (/(PA(1) - P0(1)), (PA(2) - P0(2))/)
               end if

               dist(p_i) = sqrt(v1(1)**2 + v1(2)**2)

               ! Find if the horizontal ray on right-side intersects with body
               miny = min(PA(2), PB(2))
               maxy = max(PA(2), PB(2))

               if (ycell .gt. miny .and. ycell .lt. maxy) then

                  ! Method #1 use ratios to divide the current panel using
                  ! y intersection and find x

                  mratio = PA(2) - ycell
                  nratio = ycell - PB(2)
                  xit = (mratio*PB(1) + nratio*PA(1))/(mratio + nratio)

                  ! Method #2 use the equation of line instead

                  !mratio = (PB(2)-PA(2))/(PB(1)-PA(1))
                  !xit = PA(1) + (ycell - PA(2))/mratio

                  ! Check to make sure that the intersection is on the right

                  if (xit .ge. xcell) countit = countit + 1

               end if

            end do

            ! Construct level set - if intersections are positive then the point
            ! lies outside (-), if odd then the point lies inside (+)
            mvd = sign(minval(dist(:)), 2*mod(countit, 2) - 1.)

            ! For first body explicitly satisfy level set, and then compare with
            ! existing level set for successive bodies
            lmda(i, j, k) = mvd
         end if
         !! ======================================= Classical approach ends ==============================================
      end do
   end do
   deallocate (dist)
   deallocate (distls)
   deallocate (PA, PB, Pcell, P0, v1, v2)

end subroutine ib_mapToGrid
