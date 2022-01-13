!!***
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!**

#include "constants.h"

subroutine HeatAD_getScalarProp(name, value)

  use Driver_interface, ONLY : Driver_abortFlash

  implicit none
  character(len=*), intent(in)  :: name
  real, intent(out)             :: value

  character(len=100)            :: errorMessage

  value = 0.
  write(errorMessage,*) '[HeatAD_getScalarProp] Unknown scalar: ',name
  call Driver_abortFlash(errorMessage)

end subroutine HeatAD_getScalarProp
