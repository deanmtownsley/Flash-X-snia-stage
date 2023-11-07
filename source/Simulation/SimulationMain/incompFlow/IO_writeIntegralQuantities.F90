!!****f* source/Simulation/SimulationMain/incompFlow/IO_writeIntegralQuantities
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
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   This interface can be implemented for all combinations of NDIM
!!   and geometry. Implementations can be expected to support
!!   at least 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flashx.dat. Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!
!!  ARGUMENTS
!!
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!! NOTES
!!
!!  In non-Cartesian geometries, not all integrated quantities may be directly
!!  meaningful as physical quantities. In particular, components of spatial
!!  vectors, like the components of momentum that appear under the "x-momentum",
!!  "y-momentum", and "z-momentum" headings, may not be meanigful as components
!!  of linear momentum if the corresponding x-, y-, z-coordinates are not
!!  linear.
!!
!! SEE ALSO
!!  Grid_getCellVolumes
!!***

!! DEVNOTE (10/23/2023): This is a stub to suppress source/IO/IOMain/IO_writeIntegralQuantities,
!!                       which writes a ".dat" file with integral quantities. This subroutine
!!                       maybe useful in the future at which point it should be populated based
!!                       on the implementation in source/IO/IOMain
subroutine IO_writeIntegralQuantities(isFirst, simTime)
   implicit none
   integer, intent(in) :: isFirst
   real, intent(in) :: simTime
end subroutine IO_writeIntegralQuantities
