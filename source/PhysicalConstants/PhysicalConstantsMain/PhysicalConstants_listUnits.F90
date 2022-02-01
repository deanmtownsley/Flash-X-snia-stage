!!****if* source/PhysicalConstants/PhysicalConstantsMain/PhysicalConstants_listUnits
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
!!  PhysicalConstants_listUnits
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_listUnits(integer(in) :: fileUnit)            
!!
!! DESCRIPTION
!!
!!  Writes the units of measurement to standard out
!! 
!!
!! ARGUMENTS
!!
!!     fileUnit - file number to write to
!!
!! NOTES
!!
!!***            

subroutine PhysicalConstants_listUnits (fileUnit)
  
  use PhysicalConstants_data, ONLY    : pc_typeUnit, pc_arrayUnit,        &
       &    pc_sizeUnit,  pc_SISystem, pc_nameUnitsBase, PC_NBASEUNITS

  implicit none
  
  integer, intent(in)                 :: fileUnit
  type (pc_typeUnit)                  :: unode
  integer                             :: i 

  !  List the units of measurement.
  write(fileUnit,912)
  do i=1, pc_sizeUnit
     unode = pc_arrayUnit(i)
     write (fileUnit, 906)i, trim(unode%name), unode%cgsValue, &
          &         trim(pc_nameUnitsBase(pc_SISystem,unode%baseUnit))
  enddo
  write(fileUnit,910)
  
  return        
  !------------------------------------------------------------------------        
904 format("---------------List all Units ------------",/,              &
         &              T15,"Unit",T45,"CGS Value",T70,"Base Unit")
912 format(T15,"Unit",T45,"CGS Value",T70,"Base Unit")

906 format(I3,A20,T40,1P,G15.5,T60,A20)
910 format("-----------End of Units--------------------")

end subroutine PhysicalConstants_listUnits
