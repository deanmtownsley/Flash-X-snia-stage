module sim_local_interface

  interface sim_interpolate1dWd
     subroutine sim_interpolate1dWd( radius, dr, dens, temp, xc12, xne20, xne22, fluff)
     implicit none
     real, intent(in) :: radius, dr
     real, intent(out):: dens, temp, xc12, xne20, xne22
     logical, intent(out):: fluff
     end subroutine
  end interface

  interface sim_LCGRandomIterate
    subroutine sim_LCGRandomIterate(state)
    implicit none
    integer, intent(inout) :: state
    end subroutine
  end interface

end module

