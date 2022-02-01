!!****f* source/IO/IO_writeProtons
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
!!  IO_writeProtons
!!
!! SYNOPSIS
!!
!!  call IO_writeProtons (integer (in) :: numberOfProtons,
!!                        integer (in) :: tags,
!!                        real    (in) :: points,
!!                        integer (in) :: pointCount)
!!
!! DESCRIPTION
!!
!!  Writes a collection of IO protons to the HDF5 plot file into the 'ProtonData'
!!  dataset of that plot file.
!! 
!! ARGUMENTS
!!
!!   numberOfProtons : number of protons 
!!   tags            : the proton tags
!!   points          : the proton points  
!!   pointCount      : the number of points for each proton
!!
!!***

subroutine IO_writeProtons (numberOfProtons, tags, points, pointCount)

  implicit none

  integer, intent (in) :: numberOfProtons
  integer, intent (in) :: tags       (:)
  real,    intent (in) :: points     (:,:,:)
  integer, intent (in) :: pointCount (:)

end subroutine IO_writeProtons
