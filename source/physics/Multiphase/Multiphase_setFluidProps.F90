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
#include "constants.h"
subroutine Multiphase_setFluidProps(solnData, facexData, faceyData, facezData, del,&
      logc, higc)
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM), intent(in) :: del
   integer, intent(IN), dimension(MDIM) :: logc, higc
end subroutine Multiphase_setFluidProps
