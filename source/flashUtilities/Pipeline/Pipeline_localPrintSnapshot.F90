!!****f* source/flashUtilities/Pipeline/Pipeline_localPrintSnapshot
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
!!  Pipeline_localPrintSnapshot
!!
!! SYNOPSIS
!! 
!!  call Pipeline_localPrintSnapshot (subroutine, intent (in) :: printItemDetails)
!!
!! DESCRIPTION
!!
!!  The routine will print a snapshot of the current local state of the pipeline.
!!  It prints all items in each of the buffer, sending and receiving channels in
!!  the pipeline on the local processor.
!!
!! ARGUMENTS
!!
!!  printItemDetails : the routine specifying the print statements for the items
!!
!! NOTES
!!
!!  The routine 'printItemDetails' must be defined outside the pipeline unit.
!!  It will contain the specifics of the items for which the pipeline was used.
!!
!!***

subroutine Pipeline_localPrintSnapshot (printItemDetails)

  implicit none

  interface
    subroutine printItemDetails (item, itemDescription)
      real,              intent (in) :: item (:)
      character (len=*), intent (in) :: itemDescription
    end subroutine printItemDetails
  end interface

  return
end subroutine Pipeline_localPrintSnapshot
