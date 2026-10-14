!!****f* source/physics/sourceTerms/Flame/Flame_step
!!
!! NAME
!!
!!  Flame_step
!!
!! SYNOPSIS
!!
!!  call Flame_step ( real(in) :: dt )
!!
!! DESCRIPTION
!!
!!  Performs a flame step.
!!
!!  Pricipal public function for Flame unit.
!!  Evolve flame forward for all blocks in blockList by one step
!!  of size dt.
!!  Flame speed and flame effects sub-units are called within this
!!  subroutine as they are required
!!  for ADR Applies unsplit reaction-diffusion operater to FLAM_MSCALAR
!!  May or may not deposit energy, depending on which
!!  Flame_Effects module had been included
!!
!! ARGUMENTS
!!
!!           dt - the time step.
!!
!! SEE ALSO
!!
!!  see Flame_interface.F90 for possible updates
!!
!!***

! this is a stub for when the Flame Unit is not included
!
! Dean Townsley 2008
!
subroutine Flame_step( dt )    
       
  implicit none
  real,    INTENT(in)                        :: dt

  return
end subroutine Flame_step
