!!****if* source/IO/localAPI/io_ptResetNextFile
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
!!  io_ptResetNextFile
!!
!! SYNOPSIS
!!
!!  io_ptResetNextFile(real, intent(IN)  :: savednext)
!!
!! DESCRIPTION
!!
!!  This allows you to reset the next time to write a particle file to the 
!!  value sacedNext.  This is necessary to preserve the next output time
!!  across restarts.
!!
!! ARGUMENTS
!!
!!  savednext : the externally stored next time to output a particle fil
!!
!!
!!
!!***

subroutine io_ptResetNextFile(savedNext)


  implicit none
  real, intent(IN) :: savedNext


end subroutine
