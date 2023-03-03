!!****f* source/RadTrans/RadTrans_molImplicitUpdate
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
!!      RadTrans_molImplicitUpdate
!!
!!  SYNOPSIS
!!
!!      call RadTrans_molImplicitUpdate(real, intent(in) :: t
!!                                      real, intent(in) :: dt)
!!
!!  DESCRIPTION 
!!
!!      Implicitly update evolved variables from t to t+dt
!!
!!
!!  ARGUMENTS
!!
!!      t  : Current time
!!      dt : Size of the time step to take
!!
!!***
subroutine RadTrans_molImplicitUpdate(t, dt)
    implicit none

    real, intent(in) :: t, dt

    return
end subroutine RadTrans_molImplicitUpdate
