!!****f* source/RadTrans/RadTrans_molPreEvolve
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  NAME 
!!
!!      RadTrans_molPreEvolve
!!
!!  SYNOPSIS
!!
!!      call RadTrans_molPreEvolve(real, intent(in) :: t)
!!
!!  DESCRIPTION 
!!
!!      Perform any pre-evolution work that must occur after all *_init
!!      and intBlock calls (e.g. setting evolved variables from primitives)
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!
!!***
subroutine RadTrans_molPreEvolve(t)
    implicit none

    real, intent(in) :: t

    return
end subroutine RadTrans_molPreEvolve
