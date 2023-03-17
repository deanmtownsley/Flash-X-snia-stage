!!****if* RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_restrictDgData
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
!!  RadTrans_restrictDgData
!!
!! SYNOPSIS
!!
!!  call RadTrans_restrictDgData(real(IN)    :: inData(:,:,:),
!!                               real(INOUT) :: outData(:,:,:),
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
!!   xface,yface,zface : cell face coordinates corresponding to the logical region
!!                       of coarse data in in the output array.
!!
!! AUTOGENROBODOC
!!
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2023-03-15 geometry support using face coords   - Austin Harris, Klaus Weide
!!
!!***

#include "Simulation.h"

subroutine RadTrans_restrictDgData(inData,outData, xface,yface,zface)

  Use TwoMoment_MeshRefinementModule, Only : &
     CoarsenX_TwoMoment

  Use GeometryFieldsModule, Only: &
     nGF, iGF_SqrtGm
  Use GeometryComputationModule, Only: &
     ComputeGeometryX
  Use MeshModule, Only : &
     MeshType, CreateMesh, DestroyMesh
  Use UnitsModule, Only : &
     Centimeter

  use RadTrans_data, ONLY : rt_str_geometry

  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
  real,intent(IN)    :: xface(:)
  real,intent(IN),OPTIONAL :: yface(:), zface(:)

  !-----Local variables
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: i1u, j1u, k1u
  Integer :: ii, jj, kk, icc, jcc, kcc

  Integer :: iNodeX
  Integer :: iFineX, nFineX(3)

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  Real    :: U_Crse(THORNADO_FLUID_NDOF)
  Real    :: U_Fine(THORNADO_FLUID_NDOF)

  Type(MeshType) :: MeshX_Crse(3)
  Type(MeshType) :: MeshX_Fine(3)

  real :: G_Crse(THORNADO_FLUID_NDOF,1,1,1,nGF)
  real :: G_Fine(THORNADO_FLUID_NDOF,1,1,1,nGF)

  integer :: iX_B1(3), iX_E1(3)
  real :: xL_Crse(3), xR_Crse(3)
  real :: xL_Fine(3), xR_Fine(3)

  nFineX = 1
  nFineX(1:NDIM) = refine_factor !or use Thornado-native variables here

  iX_B1 = 1
  iX_E1 = 1

  ! loop over thornado elements in parent block
           ! unk offsets of first child element in child block
  do k0 = 1, size(inData,3), THORNADO_NNODESX*nFineX(3)
     do j0 = 1, size(inData,2), THORNADO_NNODESX*nFineX(2)
        do i0 = 1, size(inData,1), THORNADO_NNODESX*nFineX(1)    ! 1, 5, 9, 13

           ! get indices indices for the coarse output data array
           ! note: no skipping by twos!
           k1 = 1 + (k0-1)/nFineX(3)*K3D
           j1 = 1 + (j0-1)/nFineX(2)*K2D
           i1 = 1 + (i0-1)/nFineX(1)     ! 1, 3, 5, 7
           !! upper indices for cells in coarse thornado element
           i1u = i1 + THORNADO_NNODESX
           j1u = j1 + THORNADO_NNODESX
           k1u = k1 + THORNADO_NNODESX
           !! extents for this coarse element
           xL_Crse = (/ xface(i1 ), yface(j1 ), zface(k1 ) /) * Centimeter
           xR_Crse = (/ xface(i1u), yface(j1u), zface(k1u) /) * Centimeter

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

                    ! Calculate sqrt(Gamma) for geometry corrections
                    G_Fine = 1.0

                    !! extents for this fine element

                    ! Going to some trouble to avoid rounding differences in
                    ! some cases where faces of the fine grid are collocated
                    ! with faces fo the coarse grid... this should work
                    ! generically, i.e., even when THORNADO_NNODESX and/or
                    ! refinement factors are different from 2... may not be
                    ! worth all this trouble... - KW
                    if (icc==1) then
                       xL_Fine(1) = xL_Crse(1)
                    else if (nFineX(1) == THORNADO_NNODESX) then
                       xL_Fine(1) = xface(i1+icc-1) * Centimeter
                    else if (mod((icc-1)*THORNADO_NNODESX,nFineX(1)) == 0) then
                       xL_Fine(1) = xface(i1 + ((icc-1)*THORNADO_NNODESX)/nFineX(1)) * Centimeter
                    else
                       xL_Fine(1) = ( xL_Crse(1)*(nFineX(1)-icc+1) + xR_Crse(1)*(icc-1) ) / nFineX(1)
                    end if
                    if (icc==nFineX(1)) then
                       xR_Fine(1) = xR_Crse(1)
                    else if (nFineX(1) == THORNADO_NNODESX) then
                       xR_Fine(1) = xface(i1+icc) * Centimeter
                    else if (mod(icc*THORNADO_NNODESX,nFineX(1)) == 0) then
                       xL_Fine(1) = xface(i1 + (icc*THORNADO_NNODESX)/nFineX(1)) * Centimeter
                    else
                       xR_Fine(1) = ( xL_Crse(1)*(nFineX(1)-icc )  + xR_Crse(1)*(icc  ) ) / nFineX(1)
                    end if

                    if (jcc==1) then
                       xL_Fine(2) = xL_Crse(2)
                    else if (nFineX(2) == THORNADO_NNODESX) then
                       xL_Fine(2) = yface(j1+jcc-1) * Centimeter
                    else
                       xL_Fine(2) = ( xL_Crse(2)*(nFineX(2)-jcc+1) + xR_Crse(2)*(jcc-1) ) / nFineX(2)
                    end if
                    if (jcc==nFineX(2)) then
                       xR_Fine(2) = xR_Crse(2)
                    else if (nFineX(2) == THORNADO_NNODESX) then
                       xR_Fine(2) = yface(j1+jcc) * Centimeter
                    else
                       xR_Fine(2) = ( xL_Crse(2)*(nFineX(2)-jcc )  + xR_Crse(2)*(jcc  ) ) / nFineX(2)
                    end if

                    if (kcc==1) then
                       xL_Fine(3) = xL_Crse(3)
                    else if (nFineX(3) == THORNADO_NNODESX) then
                       xL_Fine(3) = zface(k1+kcc-1) * Centimeter
                    else
                       xL_Fine(3) = ( xL_Crse(3)*(nFineX(3)-kcc+1) + xR_Crse(3)*(kcc-1) ) / nFineX(3)
                    end if
                    if (kcc==nFineX(3)) then
                       xR_Fine(3) = xR_Crse(3)
                    else if (nFineX(3) == THORNADO_NNODESX) then
                       xR_Fine(3) = zface(k1+kcc) * Centimeter
                    else
                       xR_Fine(3) = ( xL_Crse(3)*(nFineX(3)-kcc )  + xR_Crse(3)*(kcc  ) ) / nFineX(3)
                    end if

                    ! Create DG mesh for 1 element where xL_Fine and xR_Fine are the 
                    ! lower- and upper-extents of the PM "supercell" (child) that maps to 
                    ! the thornado element
                    call CreateMesh( MeshX_Fine(1), 1, THORNADO_NNODESX, 0, xL_Fine(1), xR_Fine(1) )
                    call CreateMesh( MeshX_Fine(2), 1, THORNADO_NNODESX, 0, xL_Fine(2), xR_Fine(2) )
                    call CreateMesh( MeshX_Fine(3), 1, THORNADO_NNODESX, 0, xL_Fine(3), xR_Fine(3) )

                    call ComputeGeometryX( iX_B1, iX_E1, iX_B1, iX_E1, G_Fine, &
                       MeshX_Option = MeshX_Fine, &
                       CoordinateSystem_Option = rt_str_geometry )

                    call DestroyMesh( MeshX_Fine(1) )
                    call DestroyMesh( MeshX_Fine(2) )
                    call DestroyMesh( MeshX_Fine(3) )

                    ! grab the data from fine grid element quadrature points
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                       jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                       ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i
                       U_Fine(iNodeX) = indata(ii,jj,kk) * G_Fine(iNodeX,1,1,1,iGF_SqrtGm)
                    end do

                    ! compute contribution to parent element from child element
                    U_Crse = U_Crse + CoarsenX_TwoMoment( iFineX, U_Fine )

                 end do
              end do
           end do

           ! Calculate sqrt(Gamma) for geometry corrections
           G_Crse = 1.0

           ! Create DG mesh for 1 element where xL_Crse and xR_Crse are the 
           ! lower- and upper-extents of the PM "supercell" (parent) that maps to 
           ! the thornado element
           call CreateMesh( MeshX_Crse(1), 1, THORNADO_NNODESX, 0, xL_Crse(1), xR_Crse(1) )
           call CreateMesh( MeshX_Crse(2), 1, THORNADO_NNODESX, 0, xL_Crse(2), xR_Crse(2) )
           call CreateMesh( MeshX_Crse(3), 1, THORNADO_NNODESX, 0, xL_Crse(3), xR_Crse(3) )

           call ComputeGeometryX( iX_B1, iX_E1, iX_B1, iX_E1, G_Crse, &
              MeshX_Option = MeshX_Crse, &
              CoordinateSystem_Option = rt_str_geometry )

           call DestroyMesh( MeshX_Crse(1) )
           call DestroyMesh( MeshX_Crse(2) )
           call DestroyMesh( MeshX_Crse(3) )

           ! store the result in parent block
           do iNodeX = 1, THORNADO_FLUID_NDOF
              kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
              jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
              ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1
              outData(ii,jj,kk) = U_Crse(iNodeX) / G_Crse(iNodeX,1,1,1,iGF_SqrtGm)
           end do

        end do
     end do
  end do

end subroutine RadTrans_restrictDgData
