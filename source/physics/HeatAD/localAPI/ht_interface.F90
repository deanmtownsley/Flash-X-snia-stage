!!****if* source/physics/HeatAD/localAPI/ht_interface
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
!!  ht_interface
!!
!! SYNOPSIS
!!  use ht_interface
!!
!! DESCRIPTION
!!***

Module ht_interface

   interface
       subroutine ht_indicators(temp,ix1,ix2,jy1,jy2,kz1,kz2,minaux,maxaux)
       implicit none
       real, dimension(:,:,:), intent(in) :: temp
       integer, intent(in) :: ix1,ix2,jy1,jy2,kz1,kz2
       real, intent(inout) :: minaux, maxaux
       end subroutine
    end interface

end Module ht_interface
