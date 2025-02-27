!!****if* source/physics/Hydro/HydroMain/unsplit/hy_setMinTimeStep
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
!! NAME
!!
!!  hy_setMinTimeStep
!!
!!
!! SYNOPSIS
!!
!!  hy_setMinTimeStep(integer(in) :: blockID,
!!                        integer(in) :: i,
!!                        integer(in) :: j,
!!                        integer(in) :: k,
!!                        real(in)    :: delta,
!!                        real(in)    :: speed)
!!
!!
!! DESCRIPTION
!!
!!  Set minimum time step and its location internally.
!!
!!
!! ARGUMENTS
!!
!!  i - local index in x direction
!!  j - local index in y direction
!!  k - local index in z direction
!!  blockID - local blockID
!!  speed - local maximum wave speed
!!
!!***



subroutine hy_setMinTimeStep(blockID,i,j,k,delta,speed)

  use Hydro_data, ONLY : hy_dtmin, hy_dtminloc, hy_dtminValid, hy_dtminCfl, hy_cfl, &
                         hy_meshMe

  implicit none

  !!$ Argument list -------------------------------------
  integer, INTENT(in) :: blockID,i,j,k
  real, INTENT(in) :: delta,speed
  !!$ ---------------------------------------------------
  real :: dt_hydro

  dt_hydro = hy_cfl*delta/abs(speed)

  if( dt_hydro < hy_dtmin) then
    hy_dtmin = dt_hydro
    hy_dtminloc(1) = i
    hy_dtminloc(2) = j
    hy_dtminloc(3) = k
    hy_dtminloc(4) = blockID
    hy_dtminloc(5) = hy_meshMe
    hy_dtminCfl    = hy_cfl
    hy_dtminValid = .TRUE.
  end if

end subroutine hy_setMinTimeStep
