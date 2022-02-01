!!****if* source/PhysicalConstants/PhysicalConstantsMain/PhysicalConstants_list
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
!!  PhysicalConstants_list
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_list(integer(in) :: fileUnit)            
!!
!! DESCRIPTION
!!
!!  Writes the physical constants to standard output
!! 
!!
!! ARGUMENTS
!!
!!     fileUnit - file number to write
!!
!! NOTES
!!
!!***            

subroutine PhysicalConstants_list(fileUnit)
  
  use PhysicalConstants_data, ONLY   : pc_sizeConstant,               &
       &       pc_arrayConstant, pc_typeConstant,                                 &
       &                          pc_SISystem, pc_nameUnitsBase, PC_NBASEUNITS
  implicit none

  integer, intent(in)                :: fileUnit
  type (pc_typeConstant)             :: cnode
  integer                            :: i, m
  

  !  List the physical constants
!  write(fileUnit,910)
  write(fileUnit,912)(pc_nameUnitsBase(pc_SISystem,i),                &
       &                           i=1,PC_NBASEUNITS)
  
  do m=1, pc_sizeConstant
     cnode = pc_arrayConstant(m)
     write(fileUnit,914)m,trim(cnode%name),cnode%cgsValue,             &
          &                           (cnode%unitExponent(i),i=1,PC_NBASEUNITS)
     
  enddo
!  write(fileUnit,916)
  return           
  !------------------------------------------------------------------------           
910 format("-----------List all CGS Constants----------------")
912 format(T5,"Constant Name",T25,"Constant Value",T36,6(6X,A3))
914 format(I3,A20,1P,G15.5,T40,6(G9.2))
916 format("-----------End Constants--------------------")

end subroutine PhysicalConstants_list

