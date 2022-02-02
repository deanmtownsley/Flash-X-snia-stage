module ut_contourSurfaceInterface
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
interface
  subroutine ut_contourSurfaceAreaBlock(nlevels, isolevels, ctrData, blkLimits, blkIndex, areas)
     implicit none
#include "Simulation.h"
#include "constants.h"

     integer,                     intent(IN)  :: nlevels
     real,dimension(nlevels),     intent(IN)  :: isolevels
     real,dimension(:,:,:),       intent(IN)  :: ctrData
     integer,dimension(HIGH,MDIM),intent(IN)  :: blkLimits
     integer,                     intent(IN)  :: blkIndex
     real,dimension(nlevels),     intent(OUT) :: areas

  end subroutine ut_contourSurfaceAreaBlock
end interface
end module ut_contourSurfaceInterface
