!!
!! Dean M. Townsley 2009
!!
!! This subroutine interpolates within a 1-dimensional WD profile
!! by averaging over a cell width.  The WD profile, stored in the
!! Simulation Unit's module data area, is averaged over an interval of
!! size "dr" at radius "radius".  Returned are the density, temperature,
!! and abundances (mass fractions).  Temperature averaging is just
!! multiplicatively mass-weighted for lack of a better thing to do.
!!
!! Averaging over an interval is used to allow a graceful mapping onto
!! the unrefined grid during initial refinement.

subroutine sim_interpolate1dWd( radius, dr, dens, temp, xhe4, xc12, xo16, xne20, fluff)

  use Simulation_data, ONLY : sim_wd_dens_tab, sim_wd_temp_tab, sim_wd_he4_tab, &
                              sim_wd_c12_tab, sim_wd_o16_tab, sim_wd_ne20_tab, &
                              sim_wd_dr_inv, sim_wd_npnts, &
                              sim_fluffBoundaryRadius, sim_densFluffInner, sim_tempFluffInner, &
                              sim_densFluffOuter, sim_tempFluffOuter
  implicit none

  real, intent(in) :: radius, dr
  real, intent(out):: dens, temp, xhe4, xc12, xo16, xne20
  logical, intent(out):: fluff

  real :: left_edge, right_edge, wd_radius
  real :: shell_left, shell_right, shellvol
  real :: vol, mass, masstemp, masshe4, massc12, masso16, massne20
  integer :: imin, imax, i

  ! want to construct a mapping that is direct if the grids are the same, and
  ! otherwise averages the cells in a reasonable way
  ! to do this we assume an evenly spaced grid for which the value
  !  of each variable is constant on the intervals

  ! calculate boundaries of averaging region in grid index coordinates
  left_edge = (radius-0.5*dr) * sim_wd_dr_inv
  right_edge = (radius+0.5*dr) * sim_wd_dr_inv
  wd_radius = sim_wd_dr_inv/sim_wd_npnts

  ! if below bounds we redefine our averaging region to cover
  ! the last point on  that end
  if (floor(right_edge) <= 0) then
     left_edge = 0.0
     right_edge = 1.0
  endif
  ! above use fluff
  if (floor(left_edge) >= sim_wd_npnts) then
     ! linear ramp out to boundary radius, then uniform
     if (radius < sim_fluffBoundaryRadius) then
        dens = exp( log(sim_densFluffInner) + (log(sim_densFluffOuter)-log(sim_densFluffInner)) * &
                                               (radius-wd_radius)/(sim_fluffBoundaryRadius-wd_radius) )
        temp = exp( log(sim_tempFluffInner) + (log(sim_tempFluffOuter)-log(sim_tempFluffInner)) * &
                                               (radius-wd_radius)/(sim_fluffBoundaryRadius-wd_radius) )
     else
        dens = sim_densFluffOuter
        temp = sim_tempFluffOuter
     endif
     xhe4 = sim_wd_he4_tab(sim_wd_npnts)
     xc12 = sim_wd_c12_tab(sim_wd_npnts)
     xo16 = sim_wd_o16_tab(sim_wd_npnts)
     xne20 = sim_wd_ne20_tab(sim_wd_npnts)
     fluff=.true.
     return
  endif
  fluff = .false.

  vol = 0.0
  mass = 0.0
  masstemp = 0.0
  masshe4  = 0.0
  massc12  = 0.0
  masso16  = 0.0
  massne20 = 0.0
  if ( right_edge > real(sim_wd_npnts) ) then
     shellvol = right_edge**3 - sim_wd_npnts**3
     
     vol = vol + shellvol
     mass = mass + sim_densFluffInner*shellvol
     masstemp = masstemp + sim_densFluffInner*sim_tempFluffInner*shellvol
     masshe4  = masshe4  + sim_densFluffInner*shellvol
  endif

  ! sum through cells in 1-d grid that this region overlaps
  imin = max(0, floor(left_edge))
  imax = min(sim_wd_npnts-1, floor(right_edge))
  do i = imin, imax
     ! average over just the portion of this cell which overlaps the averaging region
     shell_left = max(real(i), left_edge)
     shell_right = min(real(i+1), right_edge)
     ! assume 1d profile has spherecial geometry to do mass averages
     ! omit prefactor constant since it will divide out
     shellvol = shell_right**3-shell_left**3

     vol  = vol  + shellvol
     mass = mass + sim_wd_dens_tab(i+1)*shellvol
     masstemp = masstemp + sim_wd_dens_tab(i+1)*sim_wd_temp_tab(i+1)*shellvol
     masshe4  = masshe4  + sim_wd_dens_tab(i+1)*sim_wd_he4_tab(i+1)*shellvol
     massc12  = massc12  + sim_wd_dens_tab(i+1)*sim_wd_c12_tab(i+1)*shellvol
     masso16  = masso16  + sim_wd_dens_tab(i+1)*sim_wd_o16_tab(i+1)*shellvol
     massne20 = massne20 + sim_wd_dens_tab(i+1)*sim_wd_ne20_tab(i+1)*shellvol
  enddo
  dens = mass/vol
  temp = masstemp/mass
  xhe4 = masshe4/mass
  xc12 = massc12/mass
  xo16 = masso16/mass
  xne20 = massne20/mass


end subroutine
