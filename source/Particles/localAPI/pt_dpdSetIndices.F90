!!****if* source/Particles/localAPI/pt_dpdSetIndices
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
!!  pt_dpdSetIndices
!!
!! SYNOPSIS
!!
!!  call pt_dpdSetIndices(integer(OUT) :: opind(3),
!!                        integer(OUT) :: npind(3),
!!                        integer(OUT) :: ovind(3),
!!                        integer(OUT) :: nvind(3),
!!                        integer(OUT) :: intvind(3),
!!                        integer(OUT) :: ofind(3),
!!                        integer(OUT) :: nfind(3))
!!
!! DESCRIPTION
!!
!! Stub 
!!
!! ARGUMENTS
!!
!!   opind : 
!!
!!   npind : 
!!
!!   ovind : 
!!
!!   nvind : 
!!
!!   intvind : 
!!
!!   ofind : 
!!
!!   nfind : 
!!
!!
!!
!!***

subroutine pt_dpdSetIndices(opind,npind,ovind,nvind,intvind,ofind,nfind)

  implicit none
#include "Simulation.h"
#include "constants.h"
  
  integer,dimension(3),INTENT(out)::opind,npind,ovind,nvind,intvind,ofind,nfind
  
 ! Set the indices for variables at n and n+1
  ! Old positions and velocities indices
  opind(:)= 0
  ovind(:)= 0
  ofind(:)= 0
  ! N positions and velocities indices
  npind(:)= 0
  nvind(:)= 0
  nfind(:)= 0
  intvind(:)= 0
  

end subroutine pt_dpdSetIndices
