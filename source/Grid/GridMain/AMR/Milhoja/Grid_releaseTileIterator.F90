!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! This is a Milhoja-specific implementation of this routine.  Please refer
!! to the documentation in this routine's stub for general interface information.
subroutine Grid_releaseTileIterator(itor)
  use Grid_iterator, ONLY : Grid_iterator_t, &
                            destroy_iterator

  implicit none

  type(Grid_iterator_t), intent(INOUT) :: itor

  call destroy_iterator(itor)
end subroutine Grid_releaseTileIterator

