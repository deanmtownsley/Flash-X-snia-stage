!!****if* source/physics/sourceTerms/Heater/HeaterMain/annSearch/htr_annBuildTree
!!
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!***

subroutine htr_annBuildTree(heater)
   use iso_c_binding
   use ANN_mod
   use ANN_types_mod
   use Heater_type, ONLY: Heater_type_t

   implicit none

   type(Heater_Type), intent(INOUT)  :: heater

   ! ANN tree variables
   integer :: siteIndex, rowSize, colSize
   real, dimension(:, :), allocatable, target :: annDset

   ! build the ANN tree for a given heater
   heater%kdTree = c_null_ptr
   colSize = heater%dims
   rowSize = heater%numSites

   allocate (annDset(rowSize, colSize))

   do siteIndex = 1, heater%numSites
      if (colSize == 2) then
         annDset(siteIndex, :) = (/heater%xSite(siteIndex), heater%ySite(siteIndex)/)
      else if (colSize == 3) then
         annDset(siteIndex, :) = (/heater%xSite(siteIndex), heater%ySite(siteIndex), heater%zSite(siteIndex)/)
      end if
   end do

   ! build the ann tree
   ! for finding nearest neighbors to get level-set value
   call ann_buildTree(rowSize, colSize, c_loc(annDset), heater%kdTree)

   deallocate (annDset)

end subroutine htr_annBuildTree
