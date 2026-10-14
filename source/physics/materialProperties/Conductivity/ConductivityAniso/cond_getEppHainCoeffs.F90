!!****if* source/physics/materialProperties/Conductivity/ConductivityAniso/cond_getEppHainCoeffs
!!
!! NAME
!!
!!  cond_getEppHainCoeffs
!!
!! SYNOPSIS
!!
!!  call cond_getEppHainCoeffs(real(IN)  :: zbar,
!!                             real(OUT) :: g0,
!!                             real(OUT) :: g0p,
!!                             real(OUT) :: g1p,
!!                             real(OUT) :: c0p,
!!                             real(OUT) :: c1p,
!!                             real(OUT) :: c2p,
!!                             real(OUT) :: g0pp,
!!                             real(OUT) :: g1pp,
!!                             real(OUT) :: c0pp,
!!                             real(OUT) :: c1pp,
!!                             real(OUT) :: c2pp)
!!
!! DESCRIPTION
!!
!! Computes the Epperlein-Haines coefficients for anisotropic electron conductivities.
!! Values are interpolated from Table IV from Epperlein & Haines, Phys. Fluids 29, 1029, 1986.
!!
!! ARGUMENTS
!!
!!   zbar             :   average ionization [unitless]
!!   g0               :   coefficient for kappa_parallel
!!   g0p, g1p         :   coefficients for kappa_perpendicular
!!   c0p, c1p, c2p    :   coefficients for kappa_perpendicular
!!   g0pp, g1pp       :   coefficients for kappa_cross
!!   c0pp, c1pp, c2pp :   coefficients for kappa_cross
!!
!! SEE ALSO
!!
!!  Conductivity 
!!***

#include "constants.h"

  subroutine cond_getEppHainCoeffs(zbar, g0, g0p, g1p, c0p, c1p, c2p, &
                                   g0pp, g1pp, c0pp, c1pp, c2pp)
    implicit none

    real, intent(in)  :: zbar
    real, intent(out) :: g0, g0p, g1p, c0p, c1p, c2p, g0pp, g1pp, c0pp, c1pp, c2pp

    integer :: i, i1, i2, nentries
    real :: myi

    real, dimension(:), allocatable :: zbar_tab, g0_tab, &
                                       g0p_tab, g1p_tab, c0p_tab, c1p_tab, c2p_tab, &
                                       g0pp_tab, c0pp_tab, c1pp_tab, c2pp_tab

    !! Set up the table
    nentries = 15
    allocate(zbar_tab(1:nentries), g0_tab(1:nentries), &
             g0p_tab(1:nentries), g1p_tab(1:nentries), &
             c0p_tab(1:nentries), c1p_tab(1:nentries), c2p_tab(1:nentries), &
             g0pp_tab(1:nentries), &
             c0pp_tab(1:nentries), c1pp_tab(1:nentries), c2pp_tab(1:nentries))
    zbar_tab = (/ 1., 2., 3., 4., 5., 6., 7., 8., 10., 12., 14., 20., 30., 60., 1000. /)
    g0_tab = (/ 3.203, 4.931, 6.115, 6.995, 7.680, 8.231, 8.685, 9.067, &
                9.673, 10.13, 10.50, 11.23, 11.90, 12.67, 13.58 /)
    g0p_tab = (/ 6.18, 9.30, 10.2, 9.14, 8.60, 8.57, 8.84, 7.93, 7.44, &
                 7.32, 7.08, 6.79, 6.74, 6.36, 6.21 /)
    g1p_tab = (/ 4.66, 3.96, 3.72, 3.60, 3.53, 3.49, 3.49, 3.43, 3.39, &
                 3.37, 3.35, 3.32, 3.30, 3.27, 3.25 /)
    c0p_tab = (/ 1.93, 1.89, 1.66, 1.31, 1.12, 1.04, 1.02, 0.875, 0.770, &
                 0.722, 0.674, 0.605, 0.566, 0.502, 0.457 /)
    c1p_tab = (/ 2.31, 3.78, 4.76, 4.63, 4.62, 4.83, 5.19, 4.74, 4.63, &
                 4.70, 4.64, 4.65, 4.81, 4.71, 4.81 /)
    c2p_tab = (/ 5.35, 7.78, 8.88, 8.80, 8.80, 8.96, 9.24, 8.84, 8.71, &
                 8.73, 8.65, 8.60, 8.66, 8.52, 8.53 /)
    g0pp_tab = (/ 4.01, 2.46, 1.13, 0.628, 0.418, 0.319, 0.268, 0.238, &
                  0.225, 0.212, 0.202, 0.200, 0.194, 0.189, 0.186 /)
    g1pp = 2.5  !! This one is constant, so does not need a table
    c0pp_tab = (/ 0.661, 0.156, 0.0442, 0.0180, 0.00963, 0.00625, 0.00461, 0.00371, &
                  0.00300, 0.00252, 0.00221, 0.00185, 0.00156, 0.00130, 0.00108 /)
    c1pp_tab = (/ 0.931, 0.398, 0.175, 0.101, 0.0702, 0.0551, 0.0465, 0.0410, &
                  0.0354, 0.0317, 0.0291, 0.0256, 0.0228, 0.0202, 0.0180 /)
    c2pp_tab = (/ 2.50, 1.71, 1.05, 0.775, 0.646, 0.578, 0.539, 0.515, & 
                  0.497, 0.482, 0.471, 0.461, 0.450, 0.440, 0.430 /)

    !! Keep constant if zbar is greater than upper table limit
    if (zbar >= zbar_tab(nentries)) then
       g0 = g0_tab(nentries)
       g0p = g0p_tab(nentries)
       g1p = g1p_tab(nentries)
       c0p = c0p_tab(nentries)
       c1p = c1p_tab(nentries)
       c2p = c2p_tab(nentries)
       g0pp = g0pp_tab(nentries)
       c0pp = c0pp_tab(nentries)
       c1pp = c1pp_tab(nentries)
       c2pp = c2pp_tab(nentries)
       return
    end if

    !! Locate position of zbar within table (myi)
    do i=1, nentries
       if (zbar < zbar_tab(i)) then
          i2 = max(i, 2)  !! leads to extrapolation if i = 1 (ie, zbar < lower table limit)
          i1 = i2 - 1
          myi = i1 + (zbar - zbar_tab(i1))/(zbar_tab(i2)-zbar_tab(i1))
          exit
       end if
    end do

    !! Perform linear interpolation/extrapolation to get coefficients
    g0 = g0_tab(i1) + (g0_tab(i2) - g0_tab(i1))*(myi - i1)
    g0p = g0p_tab(i1) + (g0p_tab(i2) - g0p_tab(i1))*(myi - i1)
    g1p = g1p_tab(i1) + (g1p_tab(i2) - g1p_tab(i1))*(myi - i1)
    c0p = c0p_tab(i1) + (c0p_tab(i2) - c0p_tab(i1))*(myi - i1)
    c1p = c1p_tab(i1) + (c1p_tab(i2) - c1p_tab(i1))*(myi - i1)
    c2p = c2p_tab(i1) + (c2p_tab(i2) - c2p_tab(i1))*(myi - i1)
    g0pp = g0pp_tab(i1) + (g0pp_tab(i2) - g0pp_tab(i1))*(myi - i1)
    c0pp = c0pp_tab(i1) + (c0pp_tab(i2) - c0pp_tab(i1))*(myi - i1)
    c1pp = c1pp_tab(i1) + (c1pp_tab(i2) - c1pp_tab(i1))*(myi - i1)
    c2pp = c2pp_tab(i1) + (c2pp_tab(i2) - c2pp_tab(i1))*(myi - i1)

    deallocate(zbar_tab, g0_tab, g0p_tab, g1p_tab, c0p_tab, c1p_tab, &
               c2p_tab, g0pp_tab, c0pp_tab, c1pp_tab, c2pp_tab)

  end subroutine cond_getEppHainCoeffs

