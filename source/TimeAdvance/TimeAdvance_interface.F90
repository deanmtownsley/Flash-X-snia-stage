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
!! @brief Public interfaces for TimeAdvance
!!
!! @details This is the header file for any time integration implementation
!!          that defines its public interfaces.

!> @ingroup TimeAdvance
!! Interfaces to TimeAdvance public procedures

Module TimeAdvance_interface

  interface
     subroutine TimeAdvance_init()
     end subroutine TimeAdvance_init
  end interface

  interface
     subroutine TimeAdvance_finalize()
     end subroutine TimeAdvance_finalize
  end interface

  interface
     subroutine TimeAdvance(dt, dtold, time)
       real, intent(IN) :: dt, dtold, time
     end subroutine TimeAdvance
  end interface

end Module TimeAdvance_interface



