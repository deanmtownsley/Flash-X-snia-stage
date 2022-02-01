!!***
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
!!
!!**

#include "constants.h"

subroutine IncompNS_setVectorProp(name, value)

  use IncompNS_data
  use Driver_interface, only: Driver_abort

  implicit none
  character(len=*), intent(in)         :: name
  real, dimension(MDIM), intent(in)    :: value

  character(len=100)                   :: errorMessage

  select case(name)
     case("Outflow_Vel_Low","outflow_vel_low","OUTFLOW_VEL_LOW")
       ins_outflowVel(LOW,:)  = value
     case("Outflow_Vel_High","outflow_vel_high","OUTFLOW_VEL_HIGH")
       ins_outflowVel(HIGH,:) = value
     case default
       write(errorMessage,*) '[IncompNS_setVectorProp] Unknown vector: ',name
       call Driver_abort(errorMessage)
  end select  

end subroutine IncompNS_setVectorProp
