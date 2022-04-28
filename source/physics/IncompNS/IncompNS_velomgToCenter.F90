!!****f* source/physics/IncompNS/IncompNS_velomgToCenter
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
!!  IncompNS_velomg2center
!!
!! SYNOPSIS
!!
!!  call IncompNS_velomg2center(integer(in) :: blockList,
!!                              integer(in) :: blockCount)
!!
!! DESCRIPTION
!!
!!   Compute cell-centered (or volume-averaged) versions of some
!!   quantities, by interpolation or averaging from face-centered
!!   versions.
!!
!! NOTES
!!
!!  To write cell centered velocities and vorticity velx,vely,velz,omgx,omgy,omgz,
!!  add REQUIRES physics/IncompNS/IncompNSExtras to Config.
!!
!!***

subroutine IncompNS_velomgToCenter()

end subroutine IncompNS_velomgToCenter
