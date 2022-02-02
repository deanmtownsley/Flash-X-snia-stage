!!****f* source/Grid/Grid_finalize
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
!!  Grid_finalize
!!
!! SYNOPSIS
!!
!!  call Grid_finalize()
!!
!! DESCRIPTION
!!
!!  Deallocates any memory that has been allocated in the Grid Unit,
!!  and calls any finalization routines of the underlying mesh package.
!!
!!
!!***


subroutine Grid_finalize()

implicit none
  return
end subroutine Grid_finalize
