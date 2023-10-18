!! source/physics/ImBound/ImBoundMain/ib_bruteForceMap
!!
!! NAME
!!
!! ib_bruteForceMap(blockCount,blockList,dt)
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
subroutine ib_bruteForceMap2D(lmda, xcenter, ycenter, dx, dy, ix1, ix2, jy1, jy2, body)

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
   integer :: i, j, k, panelIndex
   real :: xcell, ycell, zcell, mvd

   ! For the algorithm
   real, dimension(2) :: PA, PB, Pcell, P0, v1
   real, allocatable, dimension(:) :: dist
   real :: u
   integer :: nelm = 2 ! Dimension for the points, 2 for (x,y) in 2-D
   integer :: countit
   real    :: miny, maxy, mratio, nratio, xit

   ! allocating data
   allocate (dist(body%numElems))

   k = 1
   do j = jy1, jy2
      do i = ix1, ix2
         dist = 1e+200
         countit = 0   ! Counter to check no. of intersections with the body
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

         do panelIndex = 1, body%numElems ! panelIndex is short for panel_index
            ! End points for the line segment of the IB
            ! PA is on the left and PB is on the right
            PA = body%elems(panelIndex)%pA(1:2)
            PB = body%elems(panelIndex)%pB(1:2)

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

            dist(panelIndex) = sqrt(v1(1)**2 + v1(2)**2)

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
      end do
   end do

   deallocate (dist)

end subroutine ib_bruteForceMap2D

subroutine ib_bruteForceMap3D(lmda, xcenter, ycenter, zcenter, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, body)

   ! Modules Used
   use ImBound_type, ONLY: ImBound_type_t
   implicit none

   ! Arguments
   real, dimension(:, :, :), intent(inout) :: lmda
   real, dimension(:), intent(in) :: xcenter, ycenter, zcenter
   type(ImBound_type_t), intent(in) :: body
   integer, intent(in) :: ix1, ix2, jy1, jy2, kz1, kz2
   real, intent(in) :: dx, dy, dz

   ! Internal Variables
   integer :: i, j, k, panelIndex
   real :: xcell, ycell, zcell, mvd

   ! For the algorithm
   real, dimension(3) :: PA, PB, P1, P0, PC, nrm, lx, ln
   real, dimension(3) :: PP, vec, PN
   real, dimension(3) :: vecA, vecB, vecW
   real :: dotD, da, db
   real :: tempMag, tempMag1, tempMag2, tempMag3
   real, allocatable, dimension(:) :: dist
   real :: du, dn
   integer :: nelm = 3 ! Dimension for the points, 3 for (x,y,z) in 3-D
   integer :: countit
   real    :: miny, maxy, mratio, nratio, xit

   ! allocating data
   allocate (dist(body%numElems))

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2
            dist = 1e+200
            countit = 0   ! Counter to check no. of intersections with the body
            PA = 0.0
            PB = 0.0
            PC = 0.0
            P0 = 0.0
            P1 = 0.0

            ! x and y coordinates for the current grid cell
            xcell = xcenter(i)
            ycell = ycenter(j)
            zcell = zcenter(k)

            ! Grid cell point
            P1 = (/xcell, ycell, zcell/)

            do panelIndex = 1, body%numElems ! panelIndex is short for panel_index
               ! End points for the line segment of the IB
               ! PA is on the left and PB is on the right
               P0 = body%elems(panelIndex)%center
               PA = body%elems(panelIndex)%pA
               PB = body%elems(panelIndex)%pB
               PC = body%elems(panelIndex)%pC

               ! Normal to plane
               nrm = body%elems(panelIndex)%normal

               ! Direction vector for rays
               ! lx - to find intersections in x-dir
               ! ln - to find intersections in normal-dir
               lx = (/1.0, 0.0, 0.0/)
               ln = nrm

               ! Procedure to calculate intersection for lx
               vec = PA - P1
               vecA = PB - PA
               vecB = PC - PA
               dotD = DOT_PRODUCT(vecA, vecB)**2 - DOT_PRODUCT(vecA, vecA)*DOT_PRODUCT(vecB, vecB)

               if (abs(DOT_PRODUCT(lx, nrm)) .lt. 1e-13) then
                  du = 0.0

               else
                  du = DOT_PRODUCT(vec, nrm)/DOT_PRODUCT(lx, nrm)

               end if

               PP = P1 + du*lx
               vecW = PP - PA

               da = (DOT_PRODUCT(vecA, vecB)*DOT_PRODUCT(vecW, vecB) &
                     - DOT_PRODUCT(vecB, vecB)*DOT_PRODUCT(vecW, vecA))/dotD

               db = (DOT_PRODUCT(vecA, vecB)*DOT_PRODUCT(vecW, vecA) &
                     - DOT_PRODUCT(vecA, vecA)*DOT_PRODUCT(vecW, vecB))/dotD

               if (da .ge. 0.0 .and. da .le. 1.0 .and. db .ge. 0.0 .and. (da + db) .le. 1.0 .and. du .gt. 0.0) &
                  countit = countit + 1

               ! Procedure to calculate intersection for ln
               dn = DOT_PRODUCT(vec, nrm)/DOT_PRODUCT(ln, nrm)
               PN = P1 + dn*ln
               vecW = PN - PA

               da = (DOT_PRODUCT(vecA, vecB)*DOT_PRODUCT(vecW, vecB) &
                     - DOT_PRODUCT(vecB, vecB)*DOT_PRODUCT(vecW, vecA))/dotD

               db = (DOT_PRODUCT(vecA, vecB)*DOT_PRODUCT(vecW, vecA) &
                     - DOT_PRODUCT(vecA, vecA)*DOT_PRODUCT(vecW, vecB))/dotD

               ! Use normal distance if intersection point within plane else
               ! use shortest distance to vertices or center
               if (da .ge. 0.0 .and. da .le. 1.0 .and. db .ge. 0.0 .and. (da + db) .le. 1.0) then

                  tempMag = NORM2(P1 - PN)
                  dist(panelIndex) = tempMag

               else

                  tempMag1 = NORM2(P1 - PA)
                  tempMag2 = NORM2(P1 - PB)
                  tempMag3 = NORM2(P1 - PC)
                  tempMag = NORM2(P1 - P0)

                  dist(panelIndex) = minval((/tempMag, tempMag1, tempMag2, tempMag3/))

               end if

            end do

            ! Construct level set - if intersections are positive then the point
            ! lies outside (-), if odd then the point lies inside (+)
            mvd = sign(minval(dist(:)), 2*mod(countit, 2) - 1.)

            ! For first body explicitly satisfy level set, and then compare with
            ! existing level set for successive bodies
            lmda(i, j, k) = mvd

         end do
      end do
   end do

   deallocate (dist)

end subroutine ib_bruteForceMap3D
