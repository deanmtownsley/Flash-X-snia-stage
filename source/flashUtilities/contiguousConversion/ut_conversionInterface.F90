!!****ih* source/flashUtilities/contiguousConversion/ut_conversionInterface
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
!! This is the public interface file for the array conversion  
!! subroutines.
!!***
module ut_conversionInterface

  implicit none      

  interface
     subroutine ut_convertToMemoryOffset(dims, elementCoord, arrayLBound, arrayUBound, memoryOffset)
       implicit none
       integer, intent(IN) :: dims
       integer, dimension(1:dims), intent(IN) :: elementCoord, arrayUBound, arrayLBound
       integer, intent(OUT) :: memoryOffset
     end subroutine ut_convertToMemoryOffset
  end interface

  interface
     subroutine ut_convertToArrayIndicies(dims, memoryOffset, arrayLBound, arrayUBound, elementCoord)
       implicit none
       integer, intent(IN) :: dims, memoryOffset
       integer, dimension(1:dims), intent(IN) :: arrayLBound, arrayUBound
       integer, dimension(1:dims), intent(OUT) :: elementCoord
     end subroutine ut_convertToArrayIndicies
  end interface

end module ut_conversionInterface
