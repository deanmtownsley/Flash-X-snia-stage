!!****f* source/numericalTools/MoL/localAPI/ml_advance
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
!!  NAME
!!
!!      ml_advance
!!
!!  SYNOPSIS
!!
!!      call ml_advance(real, intent(in) :: t
!!                      real, intent(in) :: dt)
!!
!!  DESCRIPTION
!!
!!      Take a timestep from t to t+dt
!!
!!  ARGUMENTS
!!
!!      t  : current time
!!      dt : size of timestep
!!
!!  TODO
!!
!!      When/if subcyling is available in Flash-X, this will extend
!!      to include a level-indicator as well
!!***
subroutine ml_advance(t, dt)
    implicit none

    real, intent(in) :: t, dt

    return
end subroutine ml_advance
