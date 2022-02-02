!!****if* source/RuntimeParameters/RuntimeParametersMain/RuntimeParameters_bcast
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
!!  RuntimeParameters_bcast
!!
!! SYNOPSIS
!!
!!  RuntimeParameters_bcast()
!!
!! DESCRIPTION
!!
!! Broadcasts parameters from the current processor to the other processors.
!! Only the master processor MASTER_PE reads in runtime parameters.
!! 
!! ARGUMENTS
!!
!!        
!!
!!
!!
!!***

subroutine RuntimeParameters_bcast()

  use RuntimeParameters_data, ONLY : rp_globalMe,parameter

  implicit none
  

  call nameValueLL_bcast(parameter, rp_globalMe)

  return

end subroutine RuntimeParameters_bcast


  

