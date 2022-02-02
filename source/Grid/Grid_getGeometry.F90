!!****f* source/Grid/Grid_getGeometry
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
!!  Grid_getGeometry
!!
!! SYNOPSIS
!!
!!  call Grid_getGeometry(integer (OUT)  :: geometry)
!!               
!!  
!! DESCRIPTION 
!!
!!  Returns the global grid geometry.
!!  valid values are (CARTESIAN, POLAR, CYLINDRICAL, SPHERICAL) defined
!!  in file "constants.h"
!!
!!
!! ARGUMENTS
!!
!!  geometry - returned value
!!
!!***

subroutine Grid_getGeometry(geometry)

  implicit none

  integer, intent(OUT) :: geometry

  geometry = 0

end subroutine Grid_getGeometry
