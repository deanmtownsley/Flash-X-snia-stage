module sim_interface

   implicit none

   interface
      subroutine sim_verifySolution(t, valid, maxError)
         implicit none
         real, intent(in) :: t
         logical, intent(out) :: valid
         real, intent(out) :: maxError
      end subroutine sim_verifySolution
   end interface

end module sim_interface
