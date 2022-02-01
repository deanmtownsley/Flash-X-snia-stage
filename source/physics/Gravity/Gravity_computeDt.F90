!!****f* source/physics/Gravity/Gravity_computeDt
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
!!  Gravity_computeDt
!!  
!! SYNOPSIS
!!
!!  Gravity_computeDt(real(:,:,:,:)  :: Uin,
!!                    real (OUT)         :: dt_grav,
!!                    integer(:)(INOUT)  :: dt_minloc(5))
!!
!! DESCRIPTION
!!
!!  Compute the timestep limiter due to the gravitational solver.
!!
!! ARGUMENTS
!!
!!  dt_grav:       Will Return the limiting timestep. Should be
!!                 set to a large value (1.D99) on input.
!!  dt_minloc(5):  An array to receive information about which
!!                 processor, block, and zone was responsible
!!                 for setting the limiting timestep.  The order
!!                 is i, j, k, b, p, where (i,j,k) = zone
!!                 indices, b = local block ID, and p = PE #.
!!                 This routine should only modify these values
!!                 if it changes dt_grav.
!!
!!***

subroutine Gravity_computeDt (Uin, dt_grav, dt_minloc)

!==============================================================================

  implicit none
  
  real,dimension(:,:,:,:) :: Uin
  
  integer, intent(INOUT) ::  dt_minloc(5)
  real,intent(OUT)       ::  dt_grav
  
  dt_grav = huge(dt_grav)
  
  return

end subroutine Gravity_computeDt
