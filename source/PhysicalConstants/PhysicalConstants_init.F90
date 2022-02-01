!!****f* source/PhysicalConstants/PhysicalConstants_init
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
!!  PhysicalConstants_init
!!
!! SYNOPSIS
!!
!!  PhysicalConstants_init()
!!
!! DESCRIPTION
!!
!! This subroutine initializes the Physical Constants databases for
!! units and constants.  Must be called in the Simulation_init
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!  pc_unitsBase   [default "CGS"] set the default system of units, either "CGS" or "MKS"
!!                 CGS:  centimeters, grams, seconds, charge=esu
!!                 MKS:  meters,  kilogramks, seconds, charge=Coloumb
!!                     both systems have temperature in Kelvin
!!
!! NOTES
!!
!!***            

subroutine PhysicalConstants_init()
  implicit none 
  return
end subroutine PhysicalConstants_init

