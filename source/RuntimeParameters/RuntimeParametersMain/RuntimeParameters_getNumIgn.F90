!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_getNumIgn
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
!!  RuntimeParameters_getNumIgn
!!
!! SYNOPSIS
!!
!!  call RuntimeParameters_getNumIgn(integer(out) :: numIgnoredParams)
!!
!! DESCRIPTION
!!
!!  Returns the number of ignored runtime parameters.
!!
!!  The 'ignored' runtime parameters are those that are found 
!!  in the flash.par but are no recognized because they were
!!  not declared in any Config file included in the setup.  
!!  (Runtime parameters MUST be declared in a Config file with the 
!!  keyword PARAMETER.  They must also be given a default
!!  value in a Config file.  The flash.par file then allows the
!!  default values to be overwritten.)  
!!
!! ARGUMENTS
!!
!!  numIgnoredParams: number of times that invalid PARAMETER lines
!!                    were found while processing flash.par.
!!
!!***

subroutine RuntimeParameters_getNumIgn (numIgnoredParams)

  use RuntimeParameters_data, ONLY : rp_numIgnoredParams

  implicit none
  integer, intent(out)                      :: numIgnoredParams

  numIgnoredParams = rp_numIgnoredParams

  return
  
end subroutine RuntimeParameters_getNumIgn
