!!****if* source/numericalTools/MoL/MoLMemory/ml_memData
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
!!    ml_memData
!!
!!  SYNOPSIS
!!    use ml_memData
!!
!!  DESCRIPTION
!!    Stores data for MoLMemory
!!
!!  NOTES
!!
!!***
module ml_memData

    implicit none

    ! Ranks are (var,i,j,k,block,dataStruct)
    real, dimension(:,:,:,:,:,:), allocatable, target, save :: scratch_data

end module ml_memData
