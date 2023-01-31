subroutine sim_verifySolution(t, valid, maxError)
   implicit none

   real, intent(in) :: t
   logical, intent(out) :: valid
   real, intent(out) :: maxError

   ! The stub should break the unit test...
   valid = .false.
   maxError = 1d40
end subroutine sim_verifySolution
