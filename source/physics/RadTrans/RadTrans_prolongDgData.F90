!!****f* RadTrans/RadTrans_prolongDgData
!! NOTICE
!!  Copyright 2023 UChicago Argonne, LLC and contributors
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
!!  RadTrans_prolongDgData
!!
!! SYNOPSIS
!!
!!  call RadTrans_prolongDgData(real(IN)    :: inData(:,:,:),
!!                              real(INOUT) :: outData(:,:,:),
!!                              integer(IN) :: skip(MDIM),
!!                              integer(IN),dimension(:)     :: xface(:),
!!                              integer(IN),dimension(:),OPTIONAL :: yface(:),
!!                              integer(IN),dimension(:),OPTIONAL :: zface(:))
!!
!! DESCRIPTION
!!
!!
!! ARGUMENTS
!!
!!   inData : real input array, may be a slice corresponding to a region of cells
!!            for one variable, taken from a larger array
!!
!!   outData : real output array, may be a slice corresponding to a region of cells
!!            for one variable from a larger array
!!
!!   skip : integer array, its values should be in the range
!!            0 ...  refine_factor*THORNADO_NNODESX - 1
!!          for the NDIM active spatial directions.
!!          For each spatial direction, it indicates by how much the first output
!!          element in that direction is offset wrt the first input element.
!!
!!   xface,yface,zface : cell face coordinates corresponding to the input array.
!!
!! AUTOGENROBODOC
!!
!! HISTORY
!!
!!  2022-09-20 Created RadTrans_prolongDgData API         - Klaus Weide
!!  2022-09-22 Added skip to the interface                - Klaus Weide
!!  2023-03-14 Added face coords as arguments             - Klaus Weide
!!***

#include "constants.h"

subroutine RadTrans_prolongDgData(inData,outData,skip, xface,yface,zface)
  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
  integer,intent(IN) :: skip(MDIM)
  real,intent(IN)    :: xface(:)
  real,intent(IN),OPTIONAL :: yface(:), zface(:)
end subroutine RadTrans_prolongDgData
