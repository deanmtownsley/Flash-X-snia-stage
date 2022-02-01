!!****f* source/IO/IO_endProtonWrite
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
!!  IO_endProtonWrite
!!
!! SYNOPSIS
!!
!!  call IO_endProtonWrite ()
!!
!! DESCRIPTION
!!
!!   This subroutine is called after all of the proton data has been written
!!   to the HDF5 plot file. It simply closes the plot file.
!!
!!***

subroutine IO_endProtonWrite ()

  implicit none

end subroutine IO_endProtonWrite
