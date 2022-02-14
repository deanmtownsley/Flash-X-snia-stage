!!****f* source/IO/IO_startRayWrite
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
!!    IO_startRayWrite
!!
!! SYNOPSIS
!!
!!    IO_startRayWrite()
!!
!! DESCRIPTION
!!
!!   This routine reopens the plot file so that laser rays can be
!!   written to it. It also creates the extendible RayData dataset in
!!   the HDF5 file by calling io_h5create_raydset.
!!
!!***
subroutine IO_startRayWrite()
  implicit none

end subroutine IO_startRayWrite
