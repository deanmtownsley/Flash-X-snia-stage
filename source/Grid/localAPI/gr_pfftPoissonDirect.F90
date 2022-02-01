!!****if* source/Grid/localAPI/gr_pfftPoissonDirect
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
!!  gr_pfftPoissonDirect
!!
!!
!! SYNOPSIS
!!
!!   gr_pfftPoissonDirect(integer(IN) :: iDirection, 
!!                      integer(IN) :: solveflag, 
!!                      integer(IN) :: inSize, 
!!                      integer(IN) :: localSize(MDIM), 
!!                      integer(IN) :: globalSize(MDIM), 
!!                      integer(IN) :: transformType(MDIM)
!!                      real(IN)    :: inArray(:)
!!                      real(OUT)   :: outArray(:))
!!
!! DESCRIPTION
!!
!!   Poisson solver routine.  This module implements the 
!!   fft based method for periodic and dirichlet problems
!!   on a uniform grid.  
!!   Isolated problems are not supported
!!
!!
!! ARGUMENTS
!!
!!   iDirection  - direction of the transform, valid values are 
!!                 PFFT_FORWARD and PFFT_INVERSE 
!!   solveflag   - Indicates the boundary conditions and therefore
!!                 the solvers to be applied.
!!                 solveflag==1 => Periodic in X, Neuman in Y and Z
!!                 solveflag==2 => Periodic in X and Y, Neuman and Z
!!                 solveflag==2 => Periodic in X, Y and Z
!!   inSize      - size of inArray and outArray
!!   localSize   - the local bounds (on myPe) of the original 3D data
!!                 to be transformed
!!   globalSize  - global size of the 3D data to be transformed
!!   transformType - The type if transform to be applied along each
!!                 - of the dimensions
!!   inArray       - single dimension array containing linearized 3D data
!!                   to be transformed
!!   outArray      - single dimension array containing transformed 
!!                   linearized 3D data
!!  
!!
!!
!!
!!***


subroutine gr_pfftPoissonDirect (iDirection, solveflag, &
     inSize, localSize, globalSize, &
     transformType, inArray, outArray)
  
#include "constants.h"  
  
  implicit none
  
  integer, intent(in)    :: iDirection, solveflag, inSize  
  integer, dimension(MDIM),intent(IN) :: localSize,globalSize,transformType
  real,  dimension(inSize),intent(IN) :: inArray
  real,  dimension(inSize),intent(OUT) :: outArray
  outArray = 0.0
  return
end subroutine gr_pfftPoissonDirect


