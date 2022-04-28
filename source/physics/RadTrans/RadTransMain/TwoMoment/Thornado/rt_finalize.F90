!!****f* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/rt_finalize
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
!!  rt_finalize
!!
!! SYNOPSIS
!!
!!  call rt_finalize ()
!!
!! DESCRIPTION
!!
!!  Cleans up Thornado
!!
!! ARGUMENTS
!!
!!***
subroutine rt_finalize ()

  use rt_data, only : rt_writeTimers
  use ThornadoInitializationModule, only : FreeThornado

  implicit none

  call FreeThornado(rt_writeTimers)

  return
end subroutine rt_finalize
