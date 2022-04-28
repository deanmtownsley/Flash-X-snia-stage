!!****if* source/physics/Eos/EosMain/WeakLib/eos_wlData
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
!!  eos_helmData
!!
!!
!! SYNOPSIS
!!
!!  use eos_helmData
!!
!! DESCRIPTION
!!
!!  General parameters for EOS WeakLib
!!
!! ARGUMENTS
!!
!!
!!***

module eos_wlData

 use wlEquationOfStateTableModule

 implicit none

 character(len=80),save :: eos_file

 type(EquationOfStateTableType), target, save :: EosNewTable
 type(EquationOfStateTableType), pointer :: eos_pointer

 integer, save :: nVariables

end module eos_wlData
