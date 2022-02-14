!!****if* source/IO/IOMain/io_checkBlockShape
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
!!  io_checkBlockShape
!!
!! SYNOPSIS
!!
!!  call io_checkBlockShape(integer,intent(IN)  :: numblocks)
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   numblocks : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

subroutine io_checkBlockShape(numBlocks)
  use IO_interface, ONLY: IO_getPrevScalar
  use Driver_interface, ONLY: Driver_abort

  implicit none

  integer,intent(IN) :: numBlocks

  integer :: fileNxb, fileNyb, fileNzb

  if (numBlocks == 0) return

  call IO_getPrevScalar("nxb", fileNxb)
  call IO_getPrevScalar("nyb", fileNyb)
  call IO_getPrevScalar("nzb", fileNzb)

  if (fileNxb .NE. NXB) call Driver_abort("NXB in checkpoint does not match this executable!")
  if (fileNyb .NE. NYB) call Driver_abort("NYB in checkpoint does not match this executable!")
  if (fileNzb .NE. NZB) call Driver_abort("NZB in checkpoint does not match this executable!")

end subroutine io_checkBlockShape
