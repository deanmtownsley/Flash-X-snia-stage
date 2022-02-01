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
subroutine reconstruct(uPlus, uMinus, data1d, flat, pLo, pHi,  ind, dx)
  use Hydro_data, ONLY : hy_limRad
  use Timers_interface, ONLY : Timers_stop, Timers_start

  implicit none

#include "Simulation.h"
#include "Spark.h"
#define NRECON HY_NUM_VARS+NSPECIES+NMASS_SCALARS

  integer, intent(IN) :: ind, pLo, pHi
  real, intent(IN) :: data1d(NRECON,pLo:pHi), flat, dx
  real, dimension(NRECON), intent(OUT) :: uPlus, uMinus

  uPlus  = data1d(:,ind)
  uMinus = data1d(:,ind)

end subroutine reconstruct
