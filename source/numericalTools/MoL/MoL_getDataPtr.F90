!!****f* source/numericalTools/MoL/MoL_getDataPtr
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
!!      MoL_getDataPtr
!!
!!  SYNOPSIS
!!
!!      call MoL_getDataPtr(class(Grid_tile_t), intent(inout) :: tileDesc
!!                          real, pointer                     :: dataPtr(:,:,:,:)
!!                          integer             intent(in)    :: dataStruct)
!!
!!  DESCRIPTION
!!
!!      Request a pointer to a MoL memory block
!!
!!  ARGUMENTS
!!
!!      tileDesc   : Meta-data and data-access for the current Grid tile/block
!!      dataPtr    : Pointer that will be set to the requested MoL memory block
!!      dataStruct : Which memory data struct, which should be one of:
!!                      - MOL_EVOLVED (current evolved variables in UNK)
!!                      - MOL_INTIAL  (inital value of evolved variables)
!!                      - MOL_RHS     (current-stage RHS storage)
!!                      - Any method/scheme specific scratch index
!!
!!***
subroutine MoL_getDataPtr(tileDesc, dataPtr, dataStruct)
    use Grid_tile, only: Grid_tile_t
    
    implicit none

    class(Grid_tile_t), intent(inout) :: tileDesc
    real, dimension(:,:,:,:), pointer :: dataPtr
    integer, intent(in) :: dataStruct

    return
end subroutine MoL_getDataPtr
