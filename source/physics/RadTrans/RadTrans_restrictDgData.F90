!!****f* RadTrans/RadTrans_restrictDgData
!!
!! NAME
!!
!!  RadTrans_restrictDgData
!!
!! SYNOPSIS
!!
!!  call RadTrans_restrictDgData(real(IN)    :: inData(:,:,:),
!!                               real(INOUT) :: outData(:,:,:))
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
!!  2022-09-20 Created RadTrans_restrictDgData API        - Klaus Weide
!!***

subroutine RadTrans_restrictDgData(inData,outData)
  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
end subroutine RadTrans_restrictDgData
