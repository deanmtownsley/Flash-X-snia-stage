!!****if* source/Grid/Grid_countParticles
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
!!  Grid_countParticles
!!
!! SYNOPSIS
!!
!!  call Grid_countParticles(integer(IN)    :: props
!!                          integer(INOUT) :: localCount,
!!                          integer(IN)    :: elementTypes,
!!                 optional,integer(IN)    :: maxCount
!!                    
!!  
!! DESCRIPTION 
!!  
!!  counts the particles per type. 
!! 
!!
!! ARGUMENTS 
!!
!!  props : number of properties of each element in the dataBuf datastructure
!!
!!  localCount : While coming in it contains the current number of elements mapped to
!!                      this processor. After all the data structure movement, the number
!!                      of local elements might change, and the new value is put back into it.
!!  elementTypes  : Count of different types of elements in the simulation
!!  maxCount : This is parameter determined at runtime, and is the maximum number of local
!!                        elements that a simulation expects to have. All the arrays that hold
!!                        particles in the Particles unit are allocated based on this number.
!!                        
!! NOTES
!!   Currently this routine is called by io_ptwriteParticles (serial/parallel), Grid_moveParticles,
!!   Particles_initData 
!!***

subroutine Grid_countParticles(props, localCount,elementTypes,maxCount)
  implicit none
  integer,intent(IN) ::  maxCount,props,elementTypes
  integer, intent(INOUT) :: localCount

  return
end subroutine Grid_countParticles