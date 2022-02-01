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

subroutine IncompNS_getScalarPropReal(name, value)

  use IncompNS_data
  use Driver_interface, only: Driver_abort

  implicit none
  character(len=*), intent(in)  :: name
  real, intent(out)             :: value

  character(len=100)             :: errorMessage

  select case(name)
     case("Reynolds_Number","reynolds_number","REYNOLDS_NUMBER")
       value = ins_invReynolds
     case("Min_Divergence","min_divergence","MIN_DIVERGENCE")
       value = ins_mindiv
     case("Max_Divergence","max_divergence","MAX_DIVERGENCE")
       value = ins_maxdiv
     case default
       value = 0.
       write(errorMessage,*) '[IncompNS_getScalarProp] Unknown scalar: ',name
       call Driver_abort(errorMessage)
  end select  

end subroutine IncompNS_getScalarPropReal

subroutine IncompNS_getScalarPropLogical(name, value)

  use IncompNS_data
  use Driver_interface, only: Driver_abort

  implicit none
  character(len=*), intent(in)  :: name
  logical, intent(out)          :: value

  character(len=100)             :: errorMessage

  select case(name)
     case("Pred_Corr_Flag","pred_corr_flag","PRED_CORR_FLAG")
       value = ins_predcorrflg 
     case default
       value = .false.
       write(errorMessage,*) '[IncompNS_getScalarProp] Unknown scalar: ',name
       call Driver_abort(errorMessage)
  end select  

end subroutine IncompNS_getScalarPropLogical
