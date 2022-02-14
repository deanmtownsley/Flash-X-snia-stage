
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
module Grid_interface

!  interface Grid_getDeltas
!     subroutine Grid_getDeltas(blockId, del)
!       integer, intent(in) :: blockId
!       real, dimension(MDIM), intent(out) :: del
!     end subroutine Grid_getDeltas
!  end interface


  contains

     subroutine Grid_getDeltas(blockId, del)
       integer, intent(in) :: blockId
       real, dimension(MDIM), intent(out) :: del

       del(:) = 1.0
     end subroutine Grid_getDeltas

end module


module Driver_interface

  contains

    subroutine Driver_abort (errorMessage)
      implicit none
      character(len=*), intent(in) :: errorMessage

    print *, errorMessage
    stop
    end subroutine Driver_abort


end module
