!!****ih* source/Grid/localAPI/gr_pfftInterfaceTypeDecl
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
!!  gr_pfftInterfaceTypeDecl
!!
!! SYNOPSIS
!!
!!  use gr_pfftInterfaceTypeDecl
!!
!! DESCRIPTION
!!
!!  Contains derived data type declarations used by gr_pfftInterface.
!!  We need to place the derived data types in a module in order that 
!!  different subroutines have access to the same type declaration.
!!
!!***

module gr_pfftInterfaceTypeDecl
  implicit none
  type PossibleGrid_t
     integer :: jProcs
     integer :: kProcs
  end type PossibleGrid_t
end module gr_pfftInterfaceTypeDecl
