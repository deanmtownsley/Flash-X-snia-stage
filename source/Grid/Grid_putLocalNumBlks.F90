!!****f* source/Grid/Grid_putLocalNumBlks
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
!!  Grid_putLocalNumBlks
!!
!! SYNOPSIS
!!
!!  call Grid_putLocalNumBlks(integer(IN) :: numBlocks)
!!  
!! DESCRIPTION 
!!  Save the number of local blocks on a processor.
!!
!!  This call is only used as part of Grid initialization when
!!  restarting from a checkpoint. The call is expected to issue from
!!  code in the IO unit that deals with reading checkpoint data.
!!
!!  Only the Paramesh4 implementation of the Grid unit requires an
!!  implementation; other Grid implementations use the stub.
!!
!! ARGUMENTS
!!  numBlocks : The number of blocks currently in use on myProcessor,
!!              provided by the caller
!!
!!***

subroutine Grid_putLocalNumBlks(numBlocks)
implicit none
  integer,intent(in) :: numBlocks

  return
end subroutine Grid_putLocalNumBlks
