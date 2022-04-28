!!****f* source/Driver/Driver_driftBlock
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
!!  Driver_driftBlock
!!
!! DESCRIPTION
!!
!!  Register the content hash of the provided block data, if it is different
!!  from the last time this block was hashed then log it to file.
!!
!! ARGUMENTS
!!
!!  src_file: source file location to log in case of changed hash
!!  src_line: source line location to log in case of changed hash
!!  blk: processor local block number
!!  ptr: block data to hash, all values will be hashed so this should not include guard cells
!!  gds: grid data struct type of data
!!
!!***
subroutine Driver_driftBlock(src_file, src_line, blk, ptr, gds)
  implicit none
  character(len=*), intent(in) :: src_file
  integer, intent(in) :: src_line
  integer, intent(in) :: blk
  real, intent(in) :: ptr(:,:,:,:)
  integer, intent(in) :: gds
end subroutine Driver_driftBlock
