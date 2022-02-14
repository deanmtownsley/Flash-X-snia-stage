!!****if* source/Particles/localAPI/pt_assignWeights
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
!! NAME
!!
!!  pt_assignWeights
!!
!! SYNOPSIS
!!  pt_assignWeights( logical(in) :: fineCoarseBdry,
!!                    real(in)  :: h(1:MDIM),
!!                    real(out) :: wx(LEFT_EDGE:RIGHT_EDGE),
!!                    real(out) :: wy(LEFT_EDGE:RIGHT_EDGE),
!!                    real(out) :: wz(LEFT_EDGE:RIGHT_EDGE))
!!
!!
!! DESCRIPTION
!!   Computes weight factors for the fraction of charge that will be assigned
!!   to neighboring mesh cells based on the Cloud-In-Cell (CIC) scheme.
!!
!! ARGUMENTS
!!   fineCoarseBdry - Indicates if the block is a fine one at a fine-coarse boundary
!!   h              -  some sort of distance of coordinate
!!   wx             -  weights in the x direction
!!   wy             -  weights in the y direction
!!   wz             -  weights in the z direction
!!
!!  
!!***

subroutine pt_assignWeights(fineCoarseBdry,h,wx,wy,wz)

#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abort

  implicit none
  logical,intent(IN) :: fineCoarseBdry
  real,dimension(MDIM), intent(IN) :: h
  real,dimension(LEFT_EDGE:RIGHT_EDGE), intent(OUT) :: wx,wy,wz

  wx(:)=0.0
  wy(:)=0.0
  wz(:)=0.0

RETURN
END SUBROUTINE pt_assignWeights

