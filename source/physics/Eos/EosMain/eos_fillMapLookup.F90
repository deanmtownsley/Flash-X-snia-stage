!!****if* source/physics/Eos/EosMain/eos_fillMapLookup
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
!! NAME
!!
!!  eos_fillMapLookup
!!
!! SYNOPSIS
!!
!!  call eos_fillMapLookup()
!!
!! DESCRIPTION
!!
!!  Alleviates the expense of the eos_variableMap function call by storing 
!!  eos map data in a module level array.
!!
!! ARGUMENTS
!!
!!***

#include "constants.h"
#include "Eos_map.h"

subroutine eos_fillMapLookup()
  use Eos_data, ONLY : eos_mapLookup
  implicit none

  !Eventually we will re-use existing constants which will enumerate 1,2,3,4,5.
!!$  integer, parameter, dimension(5) :: StructLookup = (/ &
!!$       CENTER, &
!!$       FACEX, &
!!$       FACEY, &
!!$       FACEZ, &
!!$       SCRATCH &
!!$       /)
  integer :: i, j, k, dataStruct
  integer, external :: eos_variableMap

  eachStruct: do k = 1, MAX_GRID_DATA_STRUCT
     dataStruct = k
     eachDirection: do j = EOS_IN,EOS_OUT
        eachVariable: do i = 1, EOSMAP_NUM_ROLES

           eos_mapLookup(i, j, k) = &
                eos_variableMap(dataStruct, i, (j-1)) 

        end do eachVariable
     end do eachDirection
  end do eachStruct
end subroutine eos_fillMapLookup
