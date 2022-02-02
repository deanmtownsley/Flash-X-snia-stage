!!****if* source/physics/Hydro/HydroMain/M1/hy_rk_reconstruct
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

  implicit none

#include "Simulation.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: ind, pLo, pHi
  real, intent(IN) :: data1d(NRECON,pLo:pHi), flat, dx
  real, dimension(NRECON), intent(OUT) :: uPlus, uMinus
  real, dimension(NRECON) :: del_p,del_m,delbar
  real :: dupw,dloc,delta
  integer :: g

  ! First compute delbar
  del_p = data1d(:,ind+1) - data1d(:,ind  )
  del_m = data1d(:,ind  ) - data1d(:,ind-1)
  do g=1,NRECON
     delbar(g) = minmod(del_p(g), del_m(g))
  end do
  ! Interpolate to plus side of zone
  uPlus  = data1d(:,ind) + 0.5*delbar
  ! Interpolate to minus side of zone
  uMinus = data1d(:,ind) - 0.5*delbar

end subroutine reconstruct

real function mc(a,b)
  implicit none
  real :: a,b
  mc = (sign(1.0,a)+sign(1.0,b))*min(abs(a),0.25*abs(a+b),abs(b))
end function mc

real function minmod(a,b)
  implicit none
  real :: a,b
  minmod=.5 * (sign(1.,a) + sign(1.,b))*min(abs(a),abs(b))
end function minmod
