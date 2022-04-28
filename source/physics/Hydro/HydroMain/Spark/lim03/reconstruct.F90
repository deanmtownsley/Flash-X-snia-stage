!!****if* source/physics/Hydro/HydroMain/Spark/reconstruct
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
!!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!!
!!  Unless required by applicable law or agreed to in writing, software
!!  distributed under the License is distributed on an "AS IS" BASIS,
!!  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!  See the License for the specific language governing permissions and
!!  limitations under the License.
!!
!!  NAME
!!
!!  hy_rk_reconstruct
!!
!!  SYNOPSIS
!!
!!  call hy_rk_reconstruct (  )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine reconstruct(uPlus, uMinus, data1d, flat, pLo, pHi, ind, dx)
  use Hydro_data, ONLY : hy_limRad
  use Timers_interface, ONLY : Timers_stop, Timers_start

  implicit none

#include "Simulation.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: ind, pLo, pHi
  real, intent(IN) :: data1d(NRECON,pLo:pHi), flat, dx
  real, dimension(NRECON), intent(OUT) :: uPlus, uMinus
  real, dimension(NRECON) :: del_p,del_m,lim_p,lim_m,delbar
  real :: dupw,dloc,delta, lim_f
  integer :: g
  real, parameter :: epsilon = 1.0e-16
  real, parameter :: epsinv = 1.0/epsilon
  real, parameter :: onethird=1./3.
  real, dimension(NRECON) :: eta, theta, q, phi


  lim_f = 1.0/(hy_limRad*dx)**2

  ! First compute undivided differences
  del_p = data1d(:,ind+1) - data1d(:,ind  )
  del_m = data1d(:,ind  ) - data1d(:,ind-1)

  eta = (del_m**2 + del_p**2)*lim_f

  ! Plus side limiter function
  theta = del_m/(del_p+1.0e-16)
  q = (2.+theta)*onethird
  do g=1,NRECON
    phi(g) = max(0., min( q(g), max(-0.5*theta(g), &
       minval((/2.*theta(g), q(g), 1.6/))) ))
  end do
  where(eta .le. 1.-epsilon)
     lim_p = q
  elsewhere (eta .ge. 1.+epsilon)
     lim_p = phi
  elsewhere
     lim_p =  0.5*((1.0 - (eta-1.0)*epsinv)*q &
           +          (1.0 + (eta-1.0)*epsinv)*phi)
  end where

  ! Minus side limiter function
  theta = del_p/(del_m+1.0e-16)
  q = (2.+theta)*onethird
  do g=1,NRECON
    phi(g) = max(0., min( q(g), max(-0.5*theta(g), &
       minval((/2.*theta(g), q(g), 1.6/))) ))
  end do
  where (eta .le. 1.-epsilon)
     lim_m = q
  elsewhere (eta .ge. 1.+epsilon)
     lim_m = phi
  elsewhere
     lim_m =  0.5*((1.0 - (eta-1.0)*epsinv)*q &
           +          (1.0 + (eta-1.0)*epsinv)*phi)
  end where

  ! Interpolate to plus side of zone
  uPlus  = data1d(:,ind) + 0.5*lim_p*del_p
  ! Interpolate to minus side of zone
  uMinus = data1d(:,ind) - 0.5*lim_m*del_m

  ! Now ensure positivity of density and pressure
  ! Switch to second-order MINMOD for ALL variables
  if (uPlus(HY_DENS) < 0.0 .OR. uMinus(HY_DENS) < 0.0 &
       .OR. uPlus(HY_PRES) < 0.0 .OR. uMinus(HY_PRES) < 0.0 &
       .OR. uPlus(HY_RHOE) < 0.0 .OR. uMinus(HY_RHOE) < 0.0) then
     do g=1,NRECON
        delbar(g) = minmod(del_p(g), del_m(g))
     enddo
     uPlus  = data1d(:,ind) + 0.5*delbar
     uMinus = data1d(:,ind) - 0.5*delbar
  endif
end subroutine reconstruct

real function minmod(a,b)
  implicit none
  real :: a,b
  minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
end function minmod
