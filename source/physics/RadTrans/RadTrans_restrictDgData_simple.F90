!!****f* RadTrans/RadTrans_restrictDgData_simplezzz
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
!!  RadTrans_restrictDgData_simple
!!
!! SYNOPSIS
!!
!!  call RadTrans_restrictDgData(real(IN)    :: inData(:,:,:),
!!                               real(INOUT) :: outData(:,:,:))
!!
!!  call RadTrans_restrictDgData_simple(real(IN)    :: inData(:,:,:),
!!                                      real(INOUT) :: outData(:,:,:))
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
!! NOTES
!!  The specific subroutine RadTrans_restrictDgData_simple can be invoked by the
!!  generic name RadTrans_restrictDgData if the caller uses the generic interface
!!  definition in the RadTrans_interface module.
!!
!! SEE ALSO
!!  RadTrans_restrictDgData
!!  RadTrans_prolongDgData_simple
!!  RadTrans_prolongDgData
!!
!! HISTORY
!!
!!  2022-09-20 Created RadTrans_restrictDgData API        - Klaus Weide
!!  2023-03-16 Named RadTrans_restrictDgData_simple       - Klaus Weide
!!***

subroutine RadTrans_restrictDgData_simple(inData,outData)
  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
end subroutine RadTrans_restrictDgData_simple
