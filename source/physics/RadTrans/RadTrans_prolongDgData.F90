!!****f* RadTrans/RadTrans_prolongDgData
!!
!! NAME
!!
!!  RadTrans_prolongDgData
!!
!! SYNOPSIS
!!
!!  call RadTrans_prolongDgData(real(IN)    :: inData(:,:,:),
!!                              real(INOUT) :: outData(:,:,:),
!!                              integer(IN) :: skip(MDIM))
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
!! AUTOGENROBODOC
!!
!! HISTORY
!!
!!  2022-09-20 Created RadTrans_prolongDgData API         - Klaus Weide
!!  2022-09-22 Added skip to the interface                - Klaus Weide
!!***

#include "constants.h"

subroutine RadTrans_prolongDgData(inData,outData,skip)
  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
  integer,intent(IN) :: skip(MDIM)
end subroutine RadTrans_prolongDgData
