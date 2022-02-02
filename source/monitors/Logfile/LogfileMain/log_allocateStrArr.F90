!!****if* source/monitors/Logfile/LogfileMain/log_allocateStrArr
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
!!  log_allocateStrArr
!!
!! SYNOPSIS
!!
!!  log_allocateStrArr(integer, intent(in)  :: length,
!!                     integer, intent(in)  :: dim)
!!
!! DESCRIPTION
!!  Allocates a 2 dimensional array and initializes each entry to an empty string
!!
!! ARGUMENTS
!!
!!   length : size of first dimension of string to allocate
!!
!!   dim : size of second dimension of log string to allocate
!!
!!
!!
!!***

subroutine log_allocateStrArr(length, dim)

  use Logfile_data, ONLY : log_strArr
  
implicit none
  integer, intent(in) :: length, dim
  integer :: i, j
  
  allocate (log_strArr(length, dim))
  
  do i = 1, length
     do j = 1, dim
        log_strArr(i,j) = ''
     end do
  end do

end subroutine log_allocateStrArr
