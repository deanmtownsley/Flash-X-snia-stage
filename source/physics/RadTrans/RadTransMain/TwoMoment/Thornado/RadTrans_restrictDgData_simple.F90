!!****if* RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_restrictDgData_simple
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
!!  The specific subroutine implemented here can be invoked by the generic name
!!  RadTrans_restrictDgData if the caller uses the generic interface definition
!!  in the RadTrans_interface module.
!!
!! SEE ALSO
!!  RadTrans_restrictDgData
!!
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2023-03-16 Named RadTrans_restrictDgData_simple - Klaus Weide
!!
!!***

#include "Simulation.h"

subroutine RadTrans_restrictDgData_simple(inData,outData)

  Use TwoMoment_MeshRefinementModule, Only : &
     CoarsenX_TwoMoment

  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)

  !-----Local variables
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: ii, jj, kk, icc, jcc, kcc

  Integer :: iNodeX
  Integer :: iFineX, nFineX(3)

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  Real    :: U_Crse(THORNADO_FLUID_NDOF)
  Real    :: U_Fine(THORNADO_FLUID_NDOF)

  nFineX = 1
  nFineX(1:NDIM) = refine_factor !or use Thornado-native variables here

  ! loop over thornado elements in parent block
           ! unk offsets of first child element in child block
  do k0 = 1, size(inData,3), THORNADO_NNODESX*nFineX(3)
     do j0 = 1, size(inData,2), THORNADO_NNODESX*nFineX(2)
        do i0 = 1, size(inData,1), THORNADO_NNODESX*nFineX(1)    ! 1, 5, 9, 13

           U_Crse = 0.0

           ! loop over child elements for this parent element
           iFineX = 0
           do kcc = 1, nFineX(3)
              do jcc = 1, nFineX(2)
                 do icc = 1, nFineX(1)
                    iFineX = iFineX + 1

                    ! unk offsets for child block
                    k = k0 + THORNADO_NNODESX*(kcc-1)
                    j = j0 + THORNADO_NNODESX*(jcc-1)
                    i = i0 + THORNADO_NNODESX*(icc-1) ! icc = 1 : 1, 5, 9, 13
                    ! icc = 2 : 3, 7, 11, 15

                    ! grab the data from fine grid element quadrature points
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                       jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                       ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i
                       U_Fine(iNodeX) = indata(ii,jj,kk)
                    end do

                    ! compute contribution to parent element from child element
                    U_Crse = U_Crse + CoarsenX_TwoMoment( iFineX, U_Fine )

                 end do
              end do
           end do

           ! get indices indices for the coarse output data array
           ! note: no skipping by twos!
           k1 = 1 + (k0-1)/nFineX(3)*K3D
           j1 = 1 + (j0-1)/nFineX(2)*K2D
           i1 = 1 + (i0-1)/nFineX(1)     ! 1, 3, 5, 7

           ! store the result in parent block
           do iNodeX = 1, THORNADO_FLUID_NDOF
              kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
              jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
              ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1
              outData(ii,jj,kk) = U_Crse(iNodeX)
           end do

        end do
     end do
  end do

end subroutine RadTrans_restrictDgData_simple
