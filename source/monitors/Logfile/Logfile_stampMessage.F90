!!****f* source/monitors/Logfile/Logfile_stampMessage
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
!!  Logfile_stampMessage
!!
!! SYNOPSIS
!!
!!  Logfile_stampMessage(character(len=*)(in)  :: string,
!!                       logical(in),optional :: force)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!
!!   string : message to stamp to the logfile
!!
!!   force : if force=true, the Logfile is stamped no matter what myPE calls
!!           Otherwise, only the MASTER_PE can stamp the logfile
!!
!!  NOTES
!!
!!    In general, only Driver_abort should set force=true.
!!    Otherwise, extremely slow behaviour or fatal errors may occur in multi-
!!    processor runs
!!
!!***


!stamp only a string message

subroutine Logfile_stampMessage( string,force)

 
  implicit none
  character(len=*), intent(in)           :: string
  logical, intent(in), optional          :: force 

  return
  
end subroutine Logfile_stampMessage

