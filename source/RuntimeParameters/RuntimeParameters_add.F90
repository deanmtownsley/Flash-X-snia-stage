!!****f* source/RuntimeParameters/RuntimeParameters_add
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
!!  RuntimeParameters_add
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_add(   char*(in) :: name,
!!             real/int/str/log(inout) :: value) 
!!                      
!!
!! DESCRIPTION
!!
!! This function adds a parameter to the runtime parameters database.
!! RuntimeParameters_add is an overloaded routine.  Underneath the
!! hood, RuntimeParameters_add implements RuntimeParameters_addReal,
!! RuntimeParameters_addInt, RuntimeParameters_addStr and
!! RuntimeParameters_addLog.
!!
!! In general, the user would not call this routine.  Runtime parameters
!! are added to a simulation by including the in the Config file with
!! the keyword PARAMETER
!!
!! ARGUMENTS
!!
!! name:       name of parameter
!! value:      parameter value
!!
!! NOTES
!!   
!!  Because RuntimeParameters_add is an overloaded function, a user calling
!!  the routine must USE the interface RuntimeParameters_interface. 
!!
!!
!!***


   
subroutine RuntimeParameters_addReal (name, value, rwState)
implicit none
  character(len=*), intent(in)          :: name
  real, intent(in)                      :: value
  integer,OPTIONAL,intent(in)           :: rwState
end subroutine RuntimeParameters_addReal


   
subroutine RuntimeParameters_addInt (name, value, rwState)
implicit none
  character(len=*), intent(in)          :: name
  integer, intent(in)                   :: value
  integer,OPTIONAL,intent(in)           :: rwState
end subroutine RuntimeParameters_addInt


   
subroutine RuntimeParameters_addStr (name, value, rwState)
implicit none
  character(len=*), intent(in)          :: name
  character(len=*), intent(in)          :: value
  integer,OPTIONAL,intent(in)           :: rwState
end subroutine RuntimeParameters_addStr


   
subroutine RuntimeParameters_addLog (name, value, rwState)
implicit none
  character(len=*), intent(in)          :: name
  logical, intent(in)                   :: value
  integer,OPTIONAL,intent(in)           :: rwState
end subroutine RuntimeParameters_addLog



