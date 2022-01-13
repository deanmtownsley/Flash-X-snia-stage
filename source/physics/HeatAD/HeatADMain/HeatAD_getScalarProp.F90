subroutine HeatAD_getScalarProp(name, value)
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.

#include "constants.h"

  use HeatAD_data
  use Driver_interface, only: Driver_abortFlash

  implicit none
  character(len=*), intent(in)  :: name
  real, intent(out)             :: value

  character(len=100)            :: errorMessage

  select case(name)
     case("Prandtl_Number","prandtl_number","PRANDTL_NUMBER")
       value = ht_Prandtl
     case("Bulk_Temp","bulk_temp","BULK_TEMP")
       value = ht_Tbulk
     case default
       value = 0.
       write(errorMessage,*) '[HeatAD_getScalarProp] Unknown scalar: ',name
       call Driver_abortFlash(errorMessage)
  end select  

  return

end subroutine HeatAD_getScalarProp
