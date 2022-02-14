!!****f* source/IO/IO_outputFinal
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
!!  IO_outputFinal
!!
!!
!! SYNOPSIS
!!  
!!  IO_outputFinal()
!!
!! DESCRIPTION
!!
!!  This routine is called after the code has exited the main timestep
!!  loop.  It outputs the last checkpoint,plotfile and particle plotfile.
!!
!!  If particles are not included a stub (empty) routine will be called.
!!
!!
!! ARGUMENTS
!!
!!
!!
!!
!!***

subroutine IO_outputFinal()

implicit none

end subroutine IO_outputFinal
