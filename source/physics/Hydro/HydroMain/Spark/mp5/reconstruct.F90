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
!!  reconstruct
!!
!!  SYNOPSIS
!!
!!  call reconstruct (  )
!!
!!  DESCRIPTION
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine reconstruct(Uplus, Uminus, U, flat, pLo, pHi, i, dx)

  implicit none

#include "Simulation.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: i, pLo, pHi
  real, intent(IN) :: U(NRECON,pLo:pHi), flat, dx
  real, dimension(NRECON), intent(OUT) :: Uplus, Uminus
  real, dimension(NRECON) :: delbar, del_p, del_m
  integer :: g
  real, parameter :: B1 = 1./60.
  real, parameter :: B2 = 4./3.
  real, parameter :: alpha = 4.
  real, parameter :: eps = 0. !1.e-10
  real :: d2m, d2, d2p, d2mp, d2mm, uul, umd, ulc, umin, umax, ump

  ! First construct unlimited quartic polynomials to plus and minus sides of zone
  ! Suresh & Huynh (2.1)
  Uplus  = B1*( 2.*U(:,i-2) - 13.*U(:,i-1) + 47.*U(:,i) + 27.*U(:,i+1) - 3.*U(:,i+2))
  Uminus = B1*(-3.*U(:,i-2) + 27.*U(:,i-1) + 47.*U(:,i) - 13.*U(:,i+1) + 2.*U(:,i+2))

  ! Now compute limited slope interpolants for i+1/2
  del_p = U(:,i+1) - U(:,i)
  del_m = U(:,i)   - U(:,i-1)
  ! Loop over reconstructed variables and limit as needed
  do g=1,NRECON
     delbar(g) = minmod(del_p(g), del_m(g))
     ! Limit i+1/2
     ! Compute monotonic value, SH 2.12
     ump = U(g,i) + minmod(del_p(g), alpha*del_m(g))
     if ((Uplus(g)-U(g,i))*(Uplus(g)-ump) >= eps) then
        d2m = U(g,i-2) - 2.*U(g,i-1) + U(g,i)
        d2  = U(g,i-1) - 2.*U(g,i)   + U(g,i+1)
        d2p = U(g,i)   - 2.*U(g,i+1) + U(g,i+2)
        d2mp = minmod(minmod(4*d2-d2p, 4*d2p-d2), minmod(d2, d2p))
        d2mm = minmod(minmod(4*d2m-d2, 4*d2-d2m), minmod(d2m, d2))
        uul = U(g,i) + alpha*del_m(g)
        umd = 0.5*(U(g,i) + U(g,i+1)) - 0.5*d2mp
        ulc = U(g,i) + 0.5*del_m(g) + B2*d2mm
        umin = max(min(U(g,i), U(g,i+1), umd), min(U(g,i), uul, ulc))
        umax = min(max(U(g,i), U(g,i+1), umd), max(U(g,i), uul, ulc))
        Uplus(g) = median(umin, Uplus(g), umax)
        !Uplus(g) = U(g,i) + 0.5*delbar(g)
     end if
     ! Limit i-1/2
     ump = U(g,i) - minmod(del_m(g), alpha*del_p(g))
     if ((U(g,i)-Uminus(g))*(ump-Uminus(g)) >= eps) then
        uul = U(g,i) - alpha*del_p(g)
        umd = 0.5*(U(g,i) + U(g,i-1)) - 0.5*d2mm
        ulc = U(g,i) - 0.5*del_p(g) + B2*d2mp
        umin = max(min(U(g,i), U(g,i-1), umd), min(U(g,i), uul, ulc))
        umax = min(max(U(g,i), U(g,i-1), umd), max(U(g,i), uul, ulc))
        Uminus(g) = median(umin, Uminus(g), umax)
        !Uminus(g) = U(g,i) - 0.5*delbar(g)
     end if

  end do
end subroutine reconstruct

real function minmod(a,b)
  implicit none
  real :: a,b
  minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
end function minmod

real function median(a,b,c)
  implicit none
  real :: a,b,c
  median = a + minmod(b-a,c-a)
end function median
