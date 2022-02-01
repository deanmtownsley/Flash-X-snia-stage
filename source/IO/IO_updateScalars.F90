!!****f* source/IO/IO_updateScalars
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
!!  IO_updateScalars
!!
!!
!! SYNOPSIS
!!
!!  IO_updateScalars() 
!!
!!
!!
!! DESCRIPTION
!!
!!  It calls all of the Unit output routines that allow a user to write out
!!  scalar data.  The purpose is to send each unit's scalars to a scalar list
!!  to be output together in a checkpoint file.
!!
!! 
!! ARGUMENTS
!!
!!
!! NOTES
!!
!!
!!
!! SEE ALSO
!!  IO_setScalar
!!
!!
!!***

subroutine IO_updateScalars()

implicit none

end subroutine IO_updateScalars
