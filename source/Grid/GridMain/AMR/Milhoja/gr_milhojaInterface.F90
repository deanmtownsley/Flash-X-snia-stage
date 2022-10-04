!> @file
!! @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock

!> @ingroup GridMilhoja
!!
!! @brief Milhoja-specific subroutines for use only by GridMain
module gr_milhojaInterface

    interface
        subroutine gr_checkMilhojaError(routineName, MH_errorCode)
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            character(LEN=*),     intent(IN) :: routineName
            integer(MILHOJA_INT), intent(IN) :: MH_errorCode
        end subroutine gr_checkMilhojaError

        subroutine gr_fillPhysicalBcCallback(C_lo, C_hi, &
                                             C_level, &
                                             C_startVar, C_nVars) bind(c)
            use iso_c_binding,     ONLY : C_PTR
            use milhoja_types_mod, ONLY : MILHOJA_INT
            implicit none
            type(C_PTR),          intent(IN), value :: C_lo
            type(C_PTR),          intent(IN), value :: C_hi
            integer(MILHOJA_INT), intent(IN), value :: C_level
            integer(MILHOJA_INT), intent(IN), value :: C_startVar
            integer(MILHOJA_INT), intent(IN), value :: C_nVars
        end subroutine gr_fillPhysicalBcCallback

        subroutine gr_updateDataForIo()
            implicit none
        end subroutine gr_updateDataForIo

        subroutine gr_markRefineDerefineCallback(level, tags, time, &
                                                 tagval, clearval) bind(c)
            use iso_c_binding,     ONLY : C_PTR, C_CHAR
            use milhoja_types_mod, ONLY : MILHOJA_INT, MILHOJA_REAL
            implicit none
            integer(MILHOJA_INT), intent(IN), value :: level
            type(C_PTR),          intent(IN), value :: tags 
            real(MILHOJA_REAL),   intent(IN), value :: time
            character(C_CHAR),    intent(IN), value :: tagval
            character(C_CHAR),    intent(IN), value :: clearval
        end subroutine gr_markRefineDerefineCallback

        subroutine gr_writePlotfile(step)
            implicit none
            integer, intent(IN) :: step
        end subroutine gr_writePlotfile
    end interface
 
end module gr_milhojaInterface

