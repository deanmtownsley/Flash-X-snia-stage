!!****if* source/numericalTools/MoL/MoLMain/ERK/ml_finalize
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
!!      ml_finalize
!!
!!  SYNOPSIS
!!
!!      call ml_finalize()
!!
!!  DESCRIPTION
!!
!!      Finalize a method of lines unit implementation
!!
!!  ARGUMENTS
!!
!!
!!***
subroutine ml_finalize()
   use erk_data, only: erk_A, erk_b, erk_c, erk_K

   implicit none

   if (allocated(erk_A)) deallocate (erk_A)
   if (allocated(erk_b)) deallocate (erk_b)
   if (allocated(erk_c)) deallocate (erk_c)
   if (allocated(erk_K)) deallocate (erk_K)
end subroutine ml_finalize
