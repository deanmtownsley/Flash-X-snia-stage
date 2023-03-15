!!****if* RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_prolongDgData
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
!!  call RadTrans_prolongDgData(real(IN)   ,dimension(:,:,:) :: inData(:,:,:),
!!                              real(INOUT),dimension(:,:,:) :: outData(:,:,:),
!!                              integer(IN),dimension(MDIM)  :: skip(3),
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
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2022-09-22 Added skip to the interface          - Klaus Weide
!!  2023-03-14 geometry support using face coords   - Austin Harris, Klaus Weide
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine RadTrans_prolongDgData(inData,outData,skip, xface,yface,zface)

  Use TwoMoment_MeshRefinementModule, Only : &
     RefineX_TwoMoment

  Use GeometryFieldsModule, Only: &
     nGF, iGF_SqrtGm
  Use GeometryComputationModule, Only: &
     ComputeGeometryX
  Use MeshModule, Only : &
     MeshType, CreateMesh, DestroyMesh
  Use UnitsModule, Only : &
     Centimeter

  use RadTrans_data, ONLY : rt_str_geometry, rt_geometry

  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
  integer,intent(IN) :: skip(MDIM)
  real,intent(IN)    :: xface(:)
  real,intent(IN),OPTIONAL :: yface(:), zface(:)

  !-----Local variables
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: i1u, j1u, k1u
!!$  Integer :: iu, ju, ku
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

  real, parameter :: conv_x = Centimeter

  nFineX = 1
  nFineX(1:NDIM) = refine_factor !or use Thornado-native variables here

  iX_B1 = 1
  iX_E1 = 1

  ! loop over coarse (thornado) elements in parent block
  do k1 = 1, size(inData,3), THORNADO_NNODESX
     do j1 = 1, size(inData,2), THORNADO_NNODESX
        do i1 = 1, size(inData,1), THORNADO_NNODESX    ! 1, 3, 5, 7

           ! Calculate sqrt(Gamma) for geometry corrections
           G_Crse = 1.0
           !! upper indices for cells in thornado element
           i1u = i1 + THORNADO_NNODESX
           j1u = j1 + THORNADO_NNODESX
           k1u = k1 + THORNADO_NNODESX

           !! extents for this element
           xL_Crse = (/ xface(i1 ), yface(j1 ), zface(k1 ) /) * conv_x
           xR_Crse = (/ xface(i1u), yface(j1u), zface(k1u) /) * conv_x

           call CreateMesh( MeshX_Crse(1), 1, THORNADO_NNODESX, 0, xL_Crse(1), xR_Crse(1) )
           call CreateMesh( MeshX_Crse(2), 1, THORNADO_NNODESX, 0, xL_Crse(2), xR_Crse(2) )
           call CreateMesh( MeshX_Crse(3), 1, THORNADO_NNODESX, 0, xL_Crse(3), xR_Crse(3) )

           call ComputeGeometryX( iX_B1, iX_E1, iX_B1, iX_E1, G_Crse, &
              MeshX_Option = MeshX_Crse, &
              CoordinateSystem_Option = rt_str_geometry )

           do iNodeX = 1, THORNADO_FLUID_NDOF
              kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
              jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
              ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1
              U_Crse(iNodeX) = inData(ii,jj,kk) * G_Crse(iNodeX,1,1,1,iGF_SqrtGm)
           end do

           call DestroyMesh( MeshX_Crse(1) )
           call DestroyMesh( MeshX_Crse(2) )
           call DestroyMesh( MeshX_Crse(3) )

           ! offsets of first child element in output data
           k0 = 1 + (nFineX(3) * (k1-1) - skip(3)) * K3D
           j0 = 1 + (nFineX(2) * (j1-1) - skip(2)) * K2D
           i0 = 1 +  nFineX(1) * (i1-1) - skip(1)      ! 1, 5, 9, 13
!!$           kl = 1 +  nFineX(3) * (k1-1) * K3D
!!$           jl = 1 +  nFineX(2) * (j1-1) * K2D
!!$           il = 1 +  nFineX(1) * (i1-1)     ! 1, 5, 9, 13

           ! loop over fine grid element for this parent element
           iFineX = 0
           do kcc = 1, nFineX(3)
              do jcc = 1, nFineX(2)
                 do icc = 1, nFineX(1)
                    iFineX = iFineX + 1

                    ! compute fine grid element quadrature points
                    U_Fine = RefineX_TwoMoment( iFineX, U_Crse )

                    ! calculate unk indices in child block
                    k = k0 + THORNADO_NNODESX*(kcc-1)
                    j = j0 + THORNADO_NNODESX*(jcc-1)
                    i = i0 + THORNADO_NNODESX*(icc-1) ! i1 = 1 ; i0 = 1  -> 1, 3
                    ! i1 = 3 ; i0 = 5  -> 5, 7
                    ! i1 = 5 ; i0 = 9  -> 9, 11
                    ! i1 = 7 ; i0 = 13 -> 13, 15

                    ! Calculate sqrt(Gamma) for geometry corrections
                    G_Fine = 1.0
                    !! upper indices for cells in thornado element
!!$                    iu = i + THORNADO_NNODESX - 1
!!$                    ju = j + THORNADO_NNODESX - 1
!!$                    ku = k + THORNADO_NNODESX - 1

                    !! extents for this element

                    ! Going to some trouble to avoid rounding differences in
                    ! some cases where faces of the fine grid are collocated
                    ! with faces fo the coarse grid... this should work
                    ! generically, i.e., even when THORNADO_NNODESX and/or
                    ! refinement factors are different from 2... may not be
                    ! worth all this trouble... - KW
                    if (icc==1) then
                       xL_Fine(1) = xL_Crse(1)
                    else if (nFineX(1) == THORNADO_NNODESX) then
                       xL_Fine(1) = xface(i1+icc-1) * conv_x
                    else if (mod((icc-1)*THORNADO_NNODESX,nFineX(1)) == 0) then
                       xL_Fine(1) = xface(i1 + ((icc-1)*THORNADO_NNODESX)/nFineX(1)) * conv_x
                    else
                       xL_Fine(1) = ( xL_Crse(1)*(nFineX(1)-icc+1) + xR_Crse(1)*(icc-1) ) / nFineX(1)
                    end if
                    if (icc==nFineX(1)) then
                       xR_Fine(1) = xR_Crse(1)
                    else if (nFineX(1) == THORNADO_NNODESX) then
                       xR_Fine(1) = xface(i1+icc) * conv_x
                    else if (mod(icc*THORNADO_NNODESX,nFineX(1)) == 0) then
                       xL_Fine(1) = xface(i1 + (icc*THORNADO_NNODESX)/nFineX(1)) * conv_x
                    else
                       xR_Fine(1) = ( xL_Crse(1)*(nFineX(1)-icc )  + xR_Crse(1)*(icc  ) ) / nFineX(1)
                    end if

                    if (jcc==1) then
                       xL_Fine(2) = xL_Crse(2)
                    else if (nFineX(2) == THORNADO_NNODESX) then
                       xL_Fine(2) = yface(j1+jcc-1) * conv_x
                    else
                       xL_Fine(2) = ( xL_Crse(2)*(nFineX(2)-jcc+1) + xR_Crse(2)*(jcc-1) ) / nFineX(2)
                    end if
                    if (jcc==nFineX(2)) then
                       xR_Fine(2) = xR_Crse(2)
                    else if (nFineX(2) == THORNADO_NNODESX) then
                       xR_Fine(2) = yface(j1+jcc) * conv_x
                    else
                       xR_Fine(2) = ( xL_Crse(2)*(nFineX(2)-jcc )  + xR_Crse(2)*(jcc  ) ) / nFineX(2)
                    end if

                    if (kcc==1) then
                       xL_Fine(3) = xL_Crse(3)
                    else if (nFineX(3) == THORNADO_NNODESX) then
                       xL_Fine(3) = zface(k1+kcc-1) * conv_x
                    else
                       xL_Fine(3) = ( xL_Crse(3)*(nFineX(3)-kcc+1) + xR_Crse(3)*(kcc-1) ) / nFineX(3)
                    end if
                    if (kcc==nFineX(3)) then
                       xR_Fine(3) = xR_Crse(3)
                    else if (nFineX(3) == THORNADO_NNODESX) then
                       xR_Fine(3) = zface(k1+kcc) * conv_x
                    else
                       xR_Fine(3) = ( xL_Crse(3)*(nFineX(3)-kcc )  + xR_Crse(3)*(kcc  ) ) / nFineX(3)
                    end if

                    call CreateMesh( MeshX_Fine(1), 1, THORNADO_NNODESX, 0, xL_Fine(1), xR_Fine(1) )
                    call CreateMesh( MeshX_Fine(2), 1, THORNADO_NNODESX, 0, xL_Fine(2), xR_Fine(2) )
                    call CreateMesh( MeshX_Fine(3), 1, THORNADO_NNODESX, 0, xL_Fine(3), xR_Fine(3) )

                    call ComputeGeometryX( iX_B1, iX_E1, iX_B1, iX_E1, G_Fine, &
                       MeshX_Option = MeshX_Fine, &
                       CoordinateSystem_Option = rt_str_geometry )

                    ! store the result in child block
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                       jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                       ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i

                       if (      kk.GE.lbound(outData,3) .AND. kk.LE.ubound(outData,3) &
                           .AND. jj.GE.lbound(outData,2) .AND. jj.LE.ubound(outData,2) &
                           .AND. ii.GE.lbound(outData,1) .AND. ii.LE.ubound(outData,1) ) then
                          outData(ii,jj,kk) = U_Fine(iNodeX) / G_Fine(iNodeX,1,1,1,iGF_SqrtGm)
                       end if
                    end do

                    call DestroyMesh( MeshX_Fine(1) )
                    call DestroyMesh( MeshX_Fine(2) )
                    call DestroyMesh( MeshX_Fine(3) )

                 end do
              end do
           end do

        end do
     end do
  end do

end subroutine RadTrans_prolongDgData
