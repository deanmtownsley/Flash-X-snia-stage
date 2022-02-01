!!****f* source/Driver/Driver_setupParallelEnv
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
!!  Driver_setupParallelEnv
!!
!! SYNOPSIS
!!
!!  Driver_setupParallelEnv()
!!
!! DESCRIPTION
!!
!!  Initialize the parallel message-passing environment,
!!  including generation of needed communicators for all units
!!  based upon runtime information such as number of mesh copies
!!  or whether directional communicators are needed.
!!
!!  ARGUMENTS
!!  
!!
!!
!!***

subroutine Driver_setupParallelEnv ()

  implicit none             


  return
end subroutine Driver_setupParallelEnv
