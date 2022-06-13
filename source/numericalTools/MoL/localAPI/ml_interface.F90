!!****h* source/numericalTools/MoL/localAPI/ml_interface
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
!!      ml_interface
!!
!!  SYNOPSIS
!!
!!      use ml_interface
!!
!!  DESCRIPTION
!!
!!      This is the header file for the method of lines time integration unit
!!      that defines its private interface
!!
!!***
module ml_interface

    implicit none


    !! ================================= !!
    !!  Initialization and finalization  !!
    !! ================================= !!

    interface
        subroutine ml_init
        end subroutine ml_init
    end interface

    interface
        subroutine ml_finalize
        end subroutine ml_finalize
    end interface


    !! ===================== !!
    !!  Advance a time step  !!
    !! ===================== !!

    interface
        subroutine ml_advance(t, dt)
            real, intent(in) :: t, dt
        end subroutine ml_advance
    end interface
end module ml_interface
