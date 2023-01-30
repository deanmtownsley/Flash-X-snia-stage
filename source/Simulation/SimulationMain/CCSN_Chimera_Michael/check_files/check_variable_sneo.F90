!!****if* source/Simulation/SimulationMain/CCSN_Chimera_Michale/check_varibale_sneo
#include "Simulation.h"
#include "constants.h"

!!Reorder(4): Uin, Uin2
subroutine check_variable_sneo(level, lo, hi, Uin, num, Uin2)
   use Grid_interface,   ONLY : Grid_getCellCoords
   use Driver_interface, ONLY : Driver_abort

   implicit none
   integer, intent(IN)  :: level
   integer, intent(IN)  :: num
   integer, intent(IN)  :: lo(1:MDIM)
   integer, intent(IN)  :: hi(1:MDIM)
   real, dimension(:,:,:,:), pointer :: Uin    !sneo
   real, dimension(:,:,:,:), optional :: Uin2
   integer :: i, j, k
   real, allocatable :: x(:), y(:)

   allocate(x(lo(IAXIS): hi(IAXIS)))
   allocate(y(lo(JAXIS): hi(JAXIS)))
 
   call Grid_getCellCoords(IAXIS, CENTER, level, lo, hi, x)
   call Grid_getCellCoords(JAXIS, CENTER, level, lo, hi, y)

   do       k = lo(KAXIS), hi(KAXIS)
      do    j = lo(JAXIS), hi(JAXIS)
         do i = lo(IAXIS), hi(IAXIS)
              if (Uin(EINT_VAR, i, j, k) <= 1.5e1) then
                      print *, "Check num = ", num
                      print *, "Value = ", Uin(EINT_VAR, i, j, k)
                      print *, "x, y = ", x(i), y(j)
                      if (present(Uin2)) then
                          print *, "Input data: "
                          print *, "EINT_VAR = ", Uin2(EINT_VAR,i,j,k)
                      end if
                      call Driver_abort("EINT going wrong")
              end if
         end do
      end do
   end do
   deallocate(x)
   deallocate(y)
end subroutine check_variable_sneo
