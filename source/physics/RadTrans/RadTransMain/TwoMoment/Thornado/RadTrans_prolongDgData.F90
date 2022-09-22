!!****if* RadTrans/RadTransMain/TwoMoment/Thornado/RadTrans_prolongDgData
!!
!! NAME
!!
!!  RadTrans_prolongDgData
!!
!! SYNOPSIS
!!
!!  call RadTrans_prolongDgData(real(IN)   ,dimension(:,:,:) :: inData(:,:,:),
!!                              real(INOUT),dimension(:,:,:) :: outData(:,:,:),
!!                              integer(IN),dimension(MDIM)  :: skip(3))
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
!! AUTHORS
!!
!! AUTHOR: Antigoni Georgiadou     DATE: 07/20/2021
!! AUTHOR: Austin Harris           DATE: 09/16/2022
!! MODIFIED: Klaus Weide           DATE: 09/20/2022
!!  2022-09-22 Added skip to the interface          - Klaus Weide
!!
!!***

#include "Simulation.h"
#include "constants.h"

subroutine RadTrans_prolongDgData(inData,outData,skip)

  Use TwoMoment_MeshRefinementModule, Only : &
     RefineX_TwoMoment

  implicit none
  real,intent(IN)    :: inData(:,:,:)
  real,intent(INOUT) :: outData(:,:,:)
  integer,intent(IN) :: skip(MDIM)

  !-----Local variables
  Integer :: i, j, k, i1, j1, k1, i0, j0, k0
  Integer :: ii, jj, kk, icc, jcc, kcc

  Integer :: iNodeX
  Integer :: iFineX, nFineX(3)

  Integer, Parameter :: refine_factor = 2 ! Thornado assumes this for now

  Real    :: U_Crse(THORNADO_FLUID_NDOF)
  Real    :: U_Fine(THORNADO_FLUID_NDOF)

  nFineX = 1
  nFineX(1:NDIM) = refine_factor !or use Thornado-natived variables here

  ! loop over coarse (thornado) elements in parent block
  do k1 = 1, size(inData,3), THORNADO_NNODESX
     do j1 = 1, size(inData,2), THORNADO_NNODESX
        do i1 = 1, size(inData,1), THORNADO_NNODESX    ! 1, 3, 5, 7

           do iNodeX = 1, THORNADO_FLUID_NDOF
              kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k1
              jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j1
              ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i1
              U_Crse(iNodeX) = inData(ii,jj,kk)
           end do

           ! offsets of first child element in output data
           k0 = 1 + (refine_factor*(k1-1)-skip(3))*K3D
           j0 = 1 + (refine_factor*(j1-1)-skip(2))*K2D
           i0 = 1 +  refine_factor*(i1-1)-skip(1)      ! 1, 5, 9, 13

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

                    ! store the result in child block
                    do iNodeX = 1, THORNADO_FLUID_NDOF
                       kk = mod( (iNodeX-1) / THORNADO_NNODESX**2,THORNADO_NNODESX ) + k
                       jj = mod( (iNodeX-1) / THORNADO_NNODESX   ,THORNADO_NNODESX ) + j
                       ii = mod( (iNodeX-1)                      ,THORNADO_NNODESX ) + i

                       if (      kk.GE.lbound(outData,3) .AND. kk.LE.ubound(outData,3) &
                           .AND. jj.GE.lbound(outData,2) .AND. jj.LE.ubound(outData,2) &
                           .AND. ii.GE.lbound(outData,1) .AND. ii.LE.ubound(outData,1) ) then
                          outData(ii,jj,kk) = U_Fine(iNodeX)
                       end if
                    end do

                 end do
              end do
           end do

        end do
     end do
  end do

end subroutine RadTrans_prolongDgData
