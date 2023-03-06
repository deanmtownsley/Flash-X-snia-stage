!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!!
!! @brief ml_memCopy stub

!> @ingroup MoLPrivate
!!
!! @brief Copy one MoL data structure to another
!!
!! @details
!! @anchor ml_memCopy_stub
!!
!! This procedure may be used to copy between internal MoL data structures
!! or to/from UNK
!!
!! @param dst  Destintation data structure as defined in @ref Mol.h
!! @param src  Source data structure as defined in @ref Mol.h
subroutine ml_memCopy(dst, src)
   implicit none

   integer, intent(in) :: dst, src

   return
end subroutine ml_memCopy
