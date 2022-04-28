!!****if* source/Particles/localAPI/pt_advanceMC
!!
!! NAME
!!
!!  pt_advanceMC
!!
!! SYNOPSIS
!!
!!  call pt_advanceMC(real(in)   :: dtNew,
!!                    integer(in):: iSpecies)
!!
!! DESCRIPTION
!!
!!  Time advancement routine for the particle module.
!!  
!!  LEAVE THE METHOD IN PARTICLES EMPTY. The functionality is implemented in
!!  RadTrans/RadTransMain/MonteCarlo.
!!
!! ARGUMENTS
!!
!!   dtNew -- current time increment
!!   iSpecies  -- index for the particle type in pt_infoType data structure
!!
!! SIDE EFFECTS
!!
!!  Updates the POS{X,Y,Z} and VEL{X,Y,Z} properties of particles in the particles structure.
!!  Sorts particles in the particles structure by calling Grid_sortParticles.
!!
!! NOTES
!!
!!  No special handling is done for the first call - it is assumed that particle
!!  initialization fills in initial velocity components properly.
!!
!!***

!===============================================================================

subroutine pt_advanceMC (dtNew, iSpecies)
    
  
  implicit none

  
  real, INTENT(in)  :: dtNew
  integer, INTENT(in) :: ind

  
end subroutine pt_advanceMC


