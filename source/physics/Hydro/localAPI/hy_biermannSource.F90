!!****if* source/physics/Hydro/localAPI/hy_biermannSource
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
!!  hy_biermannSource
!!
!! SYNOPSIS
!!
!!  call hy_biermannSource( integer (IN) :: blockCount,
!!                         integer (IN) :: blockList(blockCount),
!!                         real    (IN) :: dt )
!!
!! DESCRIPTION
!!
!! Implement Biermann Battery Term as a source to the magnetic field. This is a stub.
!!
!! ARGUMENTS
!!
!!  blockCount -  the number of blocks in blockList
!!  blockList  -  array holding local IDs of blocks on which to advance
!!  dt         -  timestep
!!***

Subroutine hy_biermannSource ( blockCount, blockList, dt )
  
  implicit none

  ! Arguments:
  integer, intent(IN) :: blockCount
  integer, intent(IN) :: blockList(blockCount)
  real,    intent(IN) :: dt

  return

End Subroutine hy_biermannSource
