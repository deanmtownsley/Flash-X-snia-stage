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
subroutine ib_annMap2D(lmda, xcenter, ycenter, dx, dy, ix1, ix2, jy1, jy2, body)

#include "Simulation.h"
#include "constants.h"

   ! Modules Used
   use ImBound_data, ONLY: ib_annQueries, ib_annIdx
   use ib_interface, ONLY: ib_annSearchTree
   use ImBound_type, ONLY: ImBound_type_t
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
   real, dimension(2) :: PA, PB, Pcell, P0, vecI
   real, allocatable, dimension(:) :: dist
   real :: u
   integer :: annIndex
   real :: eps = 1e-13, dotNorm, magNorm, thetaNorm

   ! define ANN parameters
   allocate (dist(ib_annQueries))

   k = 1
   do j = jy1, jy2
      do i = ix1, ix2
         dist = 1e+200
         mvd = 1e+200
         PA = 0.0
         PB = 0.0
         P0 = 0.0
         Pcell = 0.0

         vecI = 0.0

         ! x and y coordinates for the current grid cell
         xcell = xcenter(i)
         ycell = ycenter(j)
         zcell = 0.0

         ! Grid cell point
         Pcell = (/xcell, ycell/)

         ! find the  nearest neighbors to compute ls value
         ib_annIdx(:) = 0
         call Timers_start("ib_annSearchTree")
         call ib_annSearchTree(body, Pcell, ib_annQueries, ib_annIdx)
         call Timers_stop("ib_annSearchTree")

         do annIndex = 1, ib_annQueries
            panelIndex = ib_annIdx(annIndex) + 1 ! need + 1 to convert c++ index to fortran

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
               vecI = (/(Pcell(1) - P0(1)), (Pcell(2) - P0(2))/)
            else
               vecI = (/(Pcell(1) - P0(1)), (Pcell(2) - P0(2))/)
            end if

            dotNorm = DOT_PRODUCT(vecI(1:2), body%elems(panelIndex)%normal(1:2))
            magNorm = NORM2(vecI(1:2))*NORM2(body%elems(panelIndex)%normal(1:2))

            thetaNorm = acos(dotNorm/(magNorm + eps))

            if ((thetaNorm <= (acos(-1.)/2 + eps)) .or. (thetaNorm > (3*acos(-1.)/2 - eps))) then
               dist(annIndex) = -NORM2(vecI(1:2))
            else
               dist(annIndex) = NORM2(vecI(1:2))
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

end subroutine ib_annMap2D

subroutine ib_annMap3D(lmda, xcenter, ycenter, zcenter, dx, dy, dz, ix1, ix2, jy1, jy2, kz1, kz2, body)

   ! Modules Used
   use ImBound_data, ONLY: ib_annQueries, ib_annIdx
   use ImBound_type, ONLY: ImBound_type_t
   use ib_interface, ONLY: ib_annSearchTree
   use Timers_interface, ONLY: Timers_start, Timers_stop
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
   real, dimension(3) :: PA, PB, Pcell, P0, PC, ln
   real, dimension(3) :: vec, PN, vecI
   real, dimension(3) :: vecA, vecB, vecW
   real :: dotD, da, db
   real, allocatable, dimension(:) :: dist
   real :: dn
   real, dimension(4) :: tempMag
   real, dimension(4, 3) :: tempVec
   integer :: vecIndex
   integer :: annIndex
   real :: eps = 1e-13, dotNorm, magNorm, thetaNorm

   ! allocating data
   allocate (dist(ib_annQueries))

   do k = kz1, kz2
      do j = jy1, jy2
         do i = ix1, ix2
            dist = 1e+200
            mvd = 1e+200
            PA = 0.0
            PB = 0.0
            PC = 0.0
            P0 = 0.0
            Pcell = 0.0

            ! x and y coordinates for the current grid cell
            xcell = xcenter(i)
            ycell = ycenter(j)
            zcell = zcenter(k)

            ! Grid cell point
            Pcell = (/xcell, ycell, zcell/)

            ! find the  nearest neighbors to compute ls value
            ib_annIdx(:) = 0
            call Timers_start("ib_annSearchTree")
            call ib_annSearchTree(body, Pcell, ib_annQueries, ib_annIdx)
            call Timers_stop("ib_annSearchTree")

            do annIndex = 1, ib_annQueries
               panelIndex = ib_annIdx(annIndex) + 1

               ! End points for the line segment of the IB
               ! PA is on the left and PB is on the right
               P0 = body%elems(panelIndex)%center
               PA = body%elems(panelIndex)%pA
               PB = body%elems(panelIndex)%pB
               PC = body%elems(panelIndex)%pC

               ! Normal to plane
               ln = body%elems(panelIndex)%normal

               ! Procedure to calculate intersection for lx
               vec = PA - Pcell
               vecA = PB - PA
               vecB = PC - PA
               dotD = DOT_PRODUCT(vecA, vecB)**2 - DOT_PRODUCT(vecA, vecA)*DOT_PRODUCT(vecB, vecB)

               ! Procedure to calculate intersection for ln
               dn = DOT_PRODUCT(vec, ln)/DOT_PRODUCT(ln, ln)
               PN = Pcell + dn*ln
               vecW = PN - PA

               da = (DOT_PRODUCT(vecA, vecB)*DOT_PRODUCT(vecW, vecB) &
                     - DOT_PRODUCT(vecB, vecB)*DOT_PRODUCT(vecW, vecA))/dotD

               db = (DOT_PRODUCT(vecA, vecB)*DOT_PRODUCT(vecW, vecA) &
                     - DOT_PRODUCT(vecA, vecA)*DOT_PRODUCT(vecW, vecB))/dotD

               ! Use normal distance if intersection point within plane else
               ! use shortest distance to vertices or center
               if (da .ge. 0.0 .and. da .le. 1.0 .and. db .ge. 0.0 .and. (da + db) .le. 1.0) then

                  vecI = Pcell - PN

               else

                  tempVec(1, :) = Pcell - PA
                  tempVec(2, :) = Pcell - PB
                  tempVec(3, :) = Pcell - PC
                  tempVec(4, :) = Pcell - P0

                  do vecIndex = 1, 4
                     tempMag(vecIndex) = NORM2(tempVec(vecIndex, :))
                  end do

                  vecI = reshape(tempVec(minloc(tempMag), :), [3])

               end if

               dotNorm = DOT_PRODUCT(vecI, ln)
               magNorm = NORM2(vecI)*NORM2(ln)

               thetaNorm = acos(dotNorm/(magNorm + eps))

               if ((thetaNorm <= (acos(-1.)/2 + eps)) .or. (thetaNorm > (3*acos(-1.)/2 - eps))) then
                  dist(annIndex) = -NORM2(vecI)
               else
                  dist(annIndex) = NORM2(vecI)
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
   end do

   deallocate (dist)

end subroutine ib_annMap3D
