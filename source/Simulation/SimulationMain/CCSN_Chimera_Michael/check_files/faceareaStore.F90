!!****if* source/Simulation/SimulationMain/CCSN_Chimera_Michale
#include "Simulation.h"
#include "constants.h"

subroutine faceareaStore(axis, level, lo, hi, Uin)
   use Driver_interface, ONLY : Driver_abort
   use Grid_interface,   ONLY : Grid_getDeltas, &
                                Grid_getCellCoords
   use Grid_data,        ONLY : gr_geometry

   implicit none
   integer, intent(IN)  :: axis
   integer, intent(IN)  :: level
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real, dimension(:,:,:,:), pointer :: Uin    !sneo

   real    :: deltas(1:MDIM)
   integer :: loCell(1:MDIM)
   integer :: hiCell(1:MDIM)

   real, allocatable :: faceCoords(:), rf(:)
#if NDIM >= 2
   real, allocatable :: thf(:)
#endif

   real    :: area, facebase
   integer :: i, j, k

   if (      (gr_geometry /= CARTESIAN) &
       .AND. (gr_geometry /= SPHERICAL   .OR. NDIM > 2) &
       .AND. (gr_geometry /= CYLINDRICAL .OR. NDIM /= 2)) then
     !areas(:, :, :) = 0.0
     call Driver_abort("[Grid_getCellFaceAreas] Not tested yet")
   else if ((axis /= IAXIS) .AND. (axis /= JAXIS) .AND. (axis /= KAXIS)) then
     call Driver_abort("[Grid_getCellFaceAreas] Invalid axis")
   end if

   call Grid_getDeltas(level, deltas)

   select case (gr_geometry)
   case (CARTESIAN)
      associate(dx => deltas(IAXIS), &
                dy => deltas(JAXIS), &
                dz => deltas(KAXIS))
#if   NDIM == 1
         area = 1.0
#elif NDIM == 2
         if      (axis == IAXIS) then
            area = dy
         else if (axis == JAXIS) then
            area = dx
         else
            ! DEV: TODO Should this set area to 1?
            call Driver_abort("[Grid_getCellFaceAreas] Invalid axis for 2D")
         end if
#elif NDIM == 3
         if      (axis == IAXIS) then
            area = dy * dz
         else if (axis == JAXIS) then
            area = dx * dz
         else if (axis == KAXIS) then
            area = dx * dy
         end if
#endif

         do       k = lo(KAXIS), hi(KAXIS)
            do    j = lo(JAXIS), hi(JAXIS)
               do i = lo(IAXIS), hi(IAXIS)
                  if (axis == IAXIS) then
                          Uin(FAAX_VAR, i, j, k) = area
                  else if (axis == JAXIS) then
                          Uin(FAAY_VAR, i, j, k) = area
                  else if (axis == KAXIS) then
                           Uin(FAAZ_VAR, i, j, k) = area
                  endif
               end do
            end do
         end do
      end associate
   case (CYLINDRICAL)
      if      (axis == IAXIS) then
         ! Get radial coordinate of face centers
         allocate(faceCoords(lo(IAXIS):hi(IAXIS)))

         ! Convert face indices to indices of cells associated with faces 
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(IAXIS, FACES, level, loCell, hiCell, faceCoords)

         associate(dz   => deltas(JAXIS), &
                   dPhi => deltas(KAXIS), &
                   r    => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                     Uin(FAAX_VAR, i, j, k) = 2.0 * PI * ABS(r(i))
                     Uin(IVAX_VAR, i, j, k) = r(i)
#elif NDIM == 2
                     Uin(FAAX_VAR, i, j, k) = 2.0 * PI * ABS(r(i)) * dz
                     Uin(IVAX_VAR, i, j, k) = r(i)
#elif NDIM == 3
                     UIN(FAAX_VAR, i, j, k) = ABS(r(i)) * dz * dPhi
                     Uin(IVAX_VAR, i, j, k) = r(i)
#endif
                  end do
               end do
            end do
         end associate
         deallocate(faceCoords)
      else if (axis == JAXIS) then
         ! Get radial coordinate of face centers
         allocate(faceCoords(lo(IAXIS):hi(IAXIS)+1))

         ! Convert face indices to indices of cells associated with faces 
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(IAXIS, FACES, level, loCell, hiCell, faceCoords)

         associate(dPhi => deltas(KAXIS), &
                   r    => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     ! DEV: TODO These can be done more simply using
                     ! the radii of the cell centers (see cell volumes)
#if   NDIM == 1
                     Uin(FAAY_VAR, i, j, k) = PI * (r(i+1)**2 - r(i)**2)
#elif NDIM == 2
                     Uin(FAAY_VAR, i, j, k) = PI * (r(i+1)**2 - r(i)**2)
#elif NDIM == 3
                     Uin(FAAY_VAR, i, j, k) = 0.5 * (r(i+1)**2 - r(i)**2) * dPhi
#endif
                  end do
               end do
            end do
         end associate
         deallocate(faceCoords)
      else if (axis == KAXIS) then
         associate(dr   => deltas(IAXIS), &
                   dz   => deltas(JAXIS))
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
#if   NDIM == 1
                     Uin(FAAZ_VAR, i, j, k) = dr
#elif NDIM == 2
                     Uin(FAAZ_VAR, i, j, k) = dr * dz
#elif NDIM == 3
                     Uin(FAAZ_VAR, i, j, k) = dr * dz
#endif
                  end do
               end do
            end do
         end associate
      end if
   case (SPHERICAL)
      if      (axis == IAXIS) then
         ! Get coordinates of faces
         allocate(faceCoords(lo(axis):hi(axis)))

         ! Convert face indices to indices of cells associated with faces
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(axis, FACES, level, loCell, hiCell, faceCoords)
#if NDIM >= 2
         allocate(thf         (lo(JAXIS):hi(JAXIS)+1))
         call Grid_getCellCoords(JAXIS, FACES, level, lo, hi, thf )
#endif

         associate(dPhi => deltas(KAXIS), &
                   r    => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     facebase = r(i) * r(i)
#if   NDIM == 1
                     Uin(FAAX_VAR, i, j, k) = facebase * 4.0 * PI
                     Uin(IVAX_VAR, i, j, k) = r(i)
#elif NDIM == 2
                     Uin(FAAX_VAR, i, j, k) = facebase * ( cos(thf(j)) - cos(thf(j+1)) ) * 2.0 * PI
                     Uin(IVAX_VAR, i, j, k) = r(i)
                     Uin(IVAY_VAR, i, j, k) = thf(j)
#elif NDIM == 3
                     Uin(FAAX_VAR, i, j, k) = facebase * ( cos(thf(j)) - cos(thf(j+1)) ) * dPhi
                     Uin(IVAX_VAR, i, j, k) = r(i)
                     Uin(IVAY_VAR, i, j, k) = thf(j)
#endif
                  end do
               end do
            end do
         end associate
         deallocate(faceCoords)
#if NDIM >= 2
         deallocate(thf)
#endif
      else if (axis == JAXIS) then
         ! Get r-coordinates of r-faces
         allocate(rf          (lo(IAXIS):hi(IAXIS)+1))
         call Grid_getCellCoords(IAXIS, FACES, level, lo, hi, rf  )
         ! Get coordinates of faces
         allocate(faceCoords(lo(axis):hi(axis)))

         ! Convert face indices to indices of cells associated with faces 
         loCell(:) = lo(:)
         hiCell(:) = hi(:)
         hiCell(axis) = hiCell(axis) - 1
         call Grid_getCellCoords(axis, FACES, level, loCell, hiCell, faceCoords)

         associate(dPhi => deltas(KAXIS), &
                   thf2  => faceCoords)
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)
                     facebase = (rf(i)+rf(i+1))*(rf(i+1)-rf(i))*0.5
#if   NDIM == 1
                     Uin(FAAY_VAR, i, j, k) = facebase * 2.0 * PI
#elif NDIM == 2
                     Uin(FAAY_VAR, i, j, k) = facebase * sin(thf2(j)) * 2.0 * PI
                     Uin(IVAY_VAR, i, j, k) = thf2(j)
#elif NDIM == 3
                     Uin(FAAY_VAR, i, j, k) = facebase * sin(thf2(j)) * dPhi
                     Uin(IVAY_VAR, i, j, k) = thf2(j)
#endif
                  end do
               end do
            end do
         end associate
         deallocate(rf)
         deallocate(faceCoords)
      else if (axis == KAXIS) then
         ! Get r-coordinates of r-faces
         allocate(rf          (lo(IAXIS):hi(IAXIS)+1))
         call Grid_getCellCoords(IAXIS, FACES, level, lo, hi, rf  )
#if NDIM >= 2
         allocate(thf         (lo(JAXIS):hi(JAXIS)+1))
         call Grid_getCellCoords(JAXIS, FACES, level, lo, hi, thf )
#endif
         associate(dthf   => deltas(JAXIS))
            do       k = lo(KAXIS), hi(KAXIS)
               do    j = lo(JAXIS), hi(JAXIS)
                  do i = lo(IAXIS), hi(IAXIS)

                     facebase = 0.5 * (rf(i+1) + rf(i)) * (rf(i+1) - rf(i))
#if NDIM == 1
                     Uin(FAAZ_VAR, i, j, k) = facebase * PI
#elif NDIM >= 2
                     Uin(FAAZ_VAR, i, j, k) = facebase * dthf
#endif
                  end do
               end do
            end do
         end associate
         deallocate(rf)
#if NDIM >= 2
         deallocate(thf)
#endif
      end if
   end select
end subroutine faceareaStore

