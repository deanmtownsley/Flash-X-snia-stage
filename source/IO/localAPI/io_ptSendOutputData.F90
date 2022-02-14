!!****if* source/IO/localAPI/io_ptSendOutputData
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
!!  io_ptSendOutputData
!!
!! SYNOPSIS
!! 
!!  io_ptSendOutputData()
!!  
!! DESCRIPTION
!!
!!  This routine adds current values of variables belonging to the IO unit that should be
!!  checkpointed to the list of variables that are written out.
!!
!!  This routine is normally called by IO_updateScalars.
!!
!! NOTES
!!
!!  Variables that begin with "io_" like io_cellDataType and
!!  io_precision are stored in the data Fortran module for the 
!!  IO unit, IO_data.  The "io_" is meant to
!!  indicate that the variable belongs to the IO Unit.
!!
!! SEE ALSO
!!
!!
!!***

subroutine io_ptSendOutputData()
  
  implicit none
  
  return
 
end subroutine io_ptSendOutputData
