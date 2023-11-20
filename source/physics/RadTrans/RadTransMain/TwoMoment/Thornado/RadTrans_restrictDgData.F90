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

subroutine RadTrans_restrictDgData(inData,outData,lmask,xface,yface,zface)

  Use TwoMoment_MeshRefinementModule, Only : &
     CoarsenX_TwoMoment
  Use ArrayUtilitiesModule, Only : &
     CreatePackIndex

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
  real,intent(IN)    :: inData(:,:,:,:)
  real,intent(INOUT) :: outData(:,:,:,:)
  logical,intent(IN) :: lmask(:)
  real,intent(IN)    :: xface(:)
  real,intent(IN),OPTIONAL :: yface(:), zface(:)

  !-----Local variables
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: i1u, j1u, k1u
  Integer :: ii, jj, kk, icc, jcc, kcc

  Integer :: iX1, iX2, iX3
  Integer :: iX1_Fine, iX2_Fine, iX3_Fine

  Integer :: iNodeX
  Integer :: iFineX, nFineX(3), nFine

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  real, allocatable :: U_Fine(:,:,:,:,:,:)
  real, allocatable :: U_Crse(:,:,:,:,:)

  real, allocatable :: G_Crse(:,:,:,:,:)
  real, allocatable :: G_Fine(:,:,:,:,:)

  integer :: nX_Crse(3)
  integer :: nX_Fine(3)

  integer :: iX_Crse_B1(3), iX_Crse_E1(3)
  integer :: iX_Fine_B1(3), iX_Fine_E1(3)

  Type(MeshType) :: MeshX_Crse(3)
  Type(MeshType) :: MeshX_Fine(3)

  integer :: lo(3), hi(3)
  real :: xL(3), xR(3)

  real, parameter :: conv_x = Centimeter

  integer :: ivar_unk2dg(size(lmask))
  integer :: ivar_dg2unk(size(lmask))
  integer :: ivar, ivar_dg, nvar_dg

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    lmask ) &
  !$OMP MAP( alloc: ivar_unk2dg, ivar_dg2unk )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA &
  !$ACC COPYIN(     lmask ) &
  !$ACC CREATE(     ivar_unk2dg, ivar_dg2unk )
#endif

  CALL CreatePackIndex( lmask, nvar_dg, ivar_unk2dg, ivar_dg2unk )

  nFineX = 1
  nFineX(1:NDIM) = refine_factor !or use Thornado-native variables here

  nFine = product( nFineX )

  lo = [ lbound(inData,2), lbound(inData,3), lbound(inData,4) ]
  hi = [ ubound(inData,2), ubound(inData,3), ubound(inData,4) ]

  !! extents for this element
  xL = (/ xface(1)          , yface(1)          , zface(1)           /) * conv_x
  xR = (/ xface(size(xface)), yface(size(yface)), zface(size(zface)) /) * conv_x

  nX_Crse = 1
  nX_Crse(1:NDIM) = ( hi(1:NDIM) - lo(1:NDIM) + 1 ) / ( THORNADO_NNODESX * nFineX(1:NDIM) )

  iX_Crse_B1 = 1
  iX_Crse_E1 = nX_Crse

  nX_Fine = 1
  nX_Fine(1:NDIM) = nX_Crse(1:NDIM) * THORNADO_NNODESX

  iX_Fine_B1 = 1
  iX_Fine_E1 = nX_Fine

  allocate( U_Crse(THORNADO_FLUID_NDOF,nX_Crse(1),nX_Crse(2),nX_Crse(3),nvar_dg) )
  allocate( U_Fine(THORNADO_FLUID_NDOF,nFine,nX_Crse(1),nX_Crse(2),nX_Crse(3),nvar_dg) )

  allocate( G_Crse(THORNADO_FLUID_NDOF,nX_Crse(1),nX_Crse(2),nX_Crse(3),nGF) )
  allocate( G_Fine(THORNADO_FLUID_NDOF,nX_Fine(1),nX_Fine(2),nX_Fine(3),nGF) )

  do k = 1, 3
     call CreateMesh( MeshX_Crse(k), nX_Crse(k), THORNADO_NNODESX, 0, xL(k), xR(k) )
     call CreateMesh( MeshX_Fine(k), nX_Fine(k), THORNADO_NNODESX, 0, xL(k), xR(k) )
  end do

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET ENTER DATA &
  !$OMP MAP( to:    inData, outData, nFineX, &
  !$OMP             nX_Crse, iX_Crse_B1, iX_Crse_E1, &
  !$OMP             nX_Fine, iX_Fine_B1, iX_Fine_E1 ) &
  !$OMP MAP( alloc: U_Crse, G_Crse, &
  !$OMP             U_Fine, G_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC ENTER DATA &
  !$ACC COPYIN(     inData, outData, nFineX, &
  !$ACC             nX_Crse, iX_Crse_B1, iX_Crse_E1, &
  !$ACC             nX_Fine, iX_Fine_B1, iX_Fine_E1 ) &
  !$ACC CREATE(     U_Crse, G_Crse, &
  !$ACC             U_Fine, G_Fine )
#endif

  ! Calculate sqrt(Gamma) for geometry corrections
  call ComputeGeometryX( iX_Crse_B1, iX_Crse_E1, iX_Crse_B1, iX_Crse_E1, G_Crse, &
     MeshX_Option = MeshX_Crse, &
     CoordinateSystem_Option = rt_str_geometry )

  call ComputeGeometryX( iX_Fine_B1, iX_Fine_E1, iX_Fine_B1, iX_Fine_E1, G_Fine, &
     MeshX_Option = MeshX_Fine, &
     CoordinateSystem_Option = rt_str_geometry )

  ! loop over fine grid elements
#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(8) &
  !$OMP PRIVATE( i1, j1, k1, iX1_Fine, iX2_Fine, iX3_Fine, iFineX, ivar, &
  !$OMP          i0, j0, k0, i, j, k, ii, jj, kk )
#elif defined( THORNADO_OACC   )
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(8) &
  !$ACC PRIVATE( i1, j1, k1, iX1_Fine, iX2_Fine, iX3_Fine, iFineX, ivar, &
  !$ACC          i0, j0, k0, i, j, k, ii, jj, kk )
#elif defined( THORNADO_OMP    )
  !$OMP PARALLEL DO COLLAPSE(8) &
  !$OMP PRIVATE( i1, j1, k1, iX1_Fine, iX2_Fine, iX3_Fine, iFineX, ivar, &
  !$OMP          i0, j0, k0, i, j, k, ii, jj, kk )
#endif
  do ivar_dg = 1, nvar_dg

     do iX3 = 1, nX_Crse(3)
        do iX2 = 1, nX_Crse(2)
           do iX1 = 1, nX_Crse(1)

              do kcc = 1, nFineX(3)
                 do jcc = 1, nFineX(2)
                    do icc = 1, nFineX(1)

                       ! store the result in child block
                       do iNodeX = 1, THORNADO_FLUID_NDOF

                          ! which child element of the parent element we are considering
                          iFineX =     icc &
                                   + ( jcc - 1 ) * nFineX(1) & 
                                   + ( kcc - 1 ) * nFineX(1) * nFineX(2)

                          ! element indexes in child block
                          iX1_Fine = ( iX1 - 1 ) * nFineX(1) + icc
                          iX2_Fine = ( iX2 - 1 ) * nFineX(2) + jcc
                          iX3_Fine = ( iX3 - 1 ) * nFineX(3) + kcc

                          ! offsets for the parent (iX1,iX2,iX3) element in unk
                          i1 = ( iX1 - 1 ) * THORNADO_NNODESX + 1
                          j1 = ( iX2 - 1 ) * THORNADO_NNODESX + 1
                          k1 = ( iX3 - 1 ) * THORNADO_NNODESX + 1

                          ! offsets of first child element in output data
                          k0 = 1 + (nFineX(3) * (k1-1)) * K3D
                          j0 = 1 + (nFineX(2) * (j1-1)) * K2D
                          i0 = 1 +  nFineX(1) * (i1-1)

                          ! calculate unk indices in child block
                          k = k0 + THORNADO_NNODESX*(kcc-1)
                          j = j0 + THORNADO_NNODESX*(jcc-1)
                          i = i0 + THORNADO_NNODESX*(icc-1)

                          kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                          jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                          ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i

                          ivar = ivar_dg2unk(ivar_dg)
                          U_Fine(iNodeX,iFineX,iX1,iX2,iX3,ivar_dg) &
                             =   inData(ivar,ii,jj,kk) &
                               * G_Fine(iNodeX,iX1_Fine,iX2_Fine,iX3_Fine,iGF_SqrtGm)

                       end do

                    end do
                 end do
              end do

           end do
        end do
     end do

  end do

  ! compute coarse grid element reconstruction
  CALL CoarsenX_TwoMoment( nX_Crse, nvar_dg, U_Fine, U_Crse )

  ! loop over coarse (thornado) elements in parent block
#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET TEAMS DISTRIBUTE PARALLEL DO SIMD COLLAPSE(5) &
  !$OMP PRIVATE( i1, j1, k1, kk, jj, ii, ivar )
#elif defined( THORNADO_OACC   )
  !$ACC PARALLEL LOOP GANG VECTOR COLLAPSE(5) &
  !$ACC PRIVATE( i1, j1, k1, kk, jj, ii, ivar )
#elif defined( THORNADO_OMP    )
  !$OMP PARALLEL DO COLLAPSE(5) &
  !$OMP PRIVATE( i1, j1, k1, kk, jj, ii, ivar )
#endif
  do ivar_dg = 1, nvar_dg

     do iX3 = 1, nX_Crse(3)
        do iX2 = 1, nX_Crse(2)
           do iX1 = 1, nX_Crse(1)

              do iNodeX = 1, THORNADO_FLUID_NDOF

                 i1 = ( iX1 - 1 ) * THORNADO_NNODESX + 1
                 j1 = ( iX2 - 1 ) * THORNADO_NNODESX + 1
                 k1 = ( iX3 - 1 ) * THORNADO_NNODESX + 1

                 kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
                 jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
                 ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1

                 ivar = ivar_dg2unk(ivar_dg)
                 outData(ivar,ii,jj,kk) &
                    =   U_Crse(iNodeX,iX1,iX2,iX3,ivar_dg) &
                      / G_Crse(iNodeX,iX1,iX2,iX3,iGF_SqrtGm)

              end do

           end do
        end do
     end do

  end do

#if   defined( THORNADO_OMP_OL )
  !$OMP TARGET EXIT DATA &
  !$OMP MAP( from:    outData ) &
  !$OMP MAP( release: inData, nFineX, lmask, ivar_unk2dg, ivar_dg2unk, &
  !$OMP               nX_Crse, iX_Crse_B1, iX_Crse_E1, U_Crse, G_Crse, &
  !$OMP               nX_Fine, iX_Fine_B1, iX_Fine_E1, U_Fine, G_Fine )
#elif defined( THORNADO_OACC   )
  !$ACC EXIT DATA &
  !$ACC COPYOUT(      outData ) &
  !$ACC DELETE(       inData, nFineX, lmask, ivar_unk2dg, ivar_dg2unk, &
  !$ACC               nX_Crse, iX_Crse_B1, iX_Crse_E1, U_Crse, G_Crse, &
  !$ACC               nX_Fine, iX_Fine_B1, iX_Fine_E1, U_Fine, G_Fine )
#endif

  do k = 1, 3
     call DestroyMesh( MeshX_Crse(k) )
     call DestroyMesh( MeshX_Fine(k) )
  end do

  deallocate( U_Crse )
  deallocate( U_Fine )
  deallocate( G_Crse )
  deallocate( G_Fine )

end subroutine RadTrans_restrictDgData
