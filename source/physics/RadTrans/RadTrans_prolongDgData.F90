!!****f* RadTrans/RadTrans_prolongDgData
!!
!! NAME
!!
!!  RadTrans_prolongDgData
!!
!! SYNOPSIS
!!
!!  call RadTrans_prolongDgData(real(IN)    :: inData(:,:,:),
!!                              real(INOUT) :: outData(:,:,:))
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
!! AUTOGENROBODOC
!!
!! HISTORY
!!
!!  2022-09-20 Created RadTrans_prolongDgData API         - Klaus Weide
!!***

subroutine RadTrans_prolongDgData(inData,outData)
  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
end subroutine RadTrans_prolongDgData
