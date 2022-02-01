!!****if* source/physics/RadTrans/localAPI/rt_interface
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
!!   rt_interface
!!
!! SYNOPSIS
!!   use rt_interface : ONLY
!!
!!  DESCRIPTION 
!!    Interface for internal RadTrans subroutines
!!
!!***
module rt_interface
  interface
     subroutine rt_init
       implicit none
     end subroutine rt_init
  end interface
  interface
     subroutine rt_finalize
       implicit none
     end subroutine rt_finalize
  end interface
end module rt_interface
