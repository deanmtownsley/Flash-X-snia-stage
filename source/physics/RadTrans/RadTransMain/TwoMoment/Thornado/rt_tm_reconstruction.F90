!!****f* source/physics/RadTrans/RadTransMain/TwoMoment/Thornado/rt_tm_reconstruction
!!
!! NAME
!!
!!  rt_tm_reconstruction
!!
!! SYNOPSIS
!!
!!  call rt_tm_reconstruction ( real(IN),pointer :: solnData(:,:,:,:), 
!!                              integer(IN) :: nX(3),
!!                              integer(IN) :: lo(MDIM),hi(MDIM),
!!                              integer(IN) :: u_lo(MDIM),u_hi(MDIM) )
!!
!! DESCRIPTION
!!
!!  Setup thornado data structures by interpolating fluid variables from Flash
!!  cell-centered values to values at DG nodes for Thornado.
!!  Copy radiation varibles to internal Thornado data structures.
!!
!! ARGUMENTS
!! solnData - pointer to cell-centered Flash data
!! nX - number of cells
!! lo - lower index limits (excluding guard cells)
!! hi - upper index limits (excluding guard cells)
!! u_lo - lower index limits (including guard cells)
!! u_hi - upper index limits (including guard cells)
!!
!!***

!!REORDER(4): solnData

subroutine rt_tm_reconstruction(solnData,nX,lo,hi,u_lo,u_hi)

#include "Flash.h"
#include "constants.h"

  use FluidFieldsModule, ONLY : uCF, iCF_D, iCF_S1, iCF_S2, iCF_S3, iCF_E, iCF_Ne
  use GeometryFieldsModule, ONLY : uGF, iGF_h_1, iGF_h_2, iGF_h_3
  use RadiationFieldsModule, ONLY : uCR
  use SubcellReconstructionModule, ONLY : ReconstructionMatrix
  use UnitsModule, ONLY : Gram, Centimeter, Second, AtomicMassUnit

  implicit none

  real, pointer, dimension(:,:,:,:) :: solnData
  integer, dimension(3), intent(in) :: nX
  integer, dimension(MDIM), intent(in) :: lo,hi,u_lo,u_hi

  integer :: i, j, k, ic, jc, kc, ii, jj, kk, n, ioff, ivar, i_sub, n_sub

  real :: D, S1, S2, S3, E, Nel
  integer :: iNode, iNodeX, iNodeE, iCF, iCR, iS, iE

  real, parameter :: conv_dens = Gram / Centimeter**3
  real, parameter :: conv_mom  = Gram / Centimeter**2 / Second
  real, parameter :: conv_enr  = Gram / Centimeter / Second**2
  real, parameter :: conv_ne   = Gram / Centimeter**3 / AtomicMassUnit
  real, parameter :: conv_J    = Gram/Second**2/Centimeter
  real, parameter :: conv_H    = Gram/Second**3

  ! Interpolate fluid vars from Flash arrays into Thornado arrays
#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(4)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG COLLAPSE(4) &
  !$ACC PRIVATE( i, j, k, D, S1, S2, S3, E, Nel ) &
  !$ACC PRESENT( ReconstructionMatrix, solnData, uCF, uGF, u_lo, u_hi, lo, hi )
#elif defined(THORNADO_OMP)
  !$OMP PARALLEL DO COLLAPSE(4)
#endif
  do kc = u_lo(KAXIS), u_hi(KAXIS)
     do jc = u_lo(JAXIS), u_hi(JAXIS)
        do ic = u_lo(IAXIS), u_hi(IAXIS)
           do iNodeX = 1, THORNADO_FLUID_NDOF

              i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
              j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
              k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

              D   = 0.0
              S1  = 0.0
              S2  = 0.0
              S3  = 0.0
              E   = 0.0
              Nel = 0.0

#if defined(THORNADO_OMP_OL)
              !$OMP PARALLEL DO SIMD &
              !$OMP PRIVATE( ii, jj, kk ) &
              !$OMP REDUCTION( +: D, S1, S2, S3, E, Nel )
#elif defined(THORNADO_OACC)
              !$ACC LOOP VECTOR &
              !$ACC PRIVATE( ii, jj, kk ) &
              !$ACC REDUCTION( +: D, S1, S2, S3, E, Nel )
#endif
              do n = 1, THORNADO_FLUID_NDOF
                 ii = mod((n-1)                    ,THORNADO_NNODESX) + i
                 jj = mod((n-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                 kk = mod((n-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k
                 D   = D   + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)
                 S1  = S1  + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(VELX_VAR,  ii,jj,kk)
                 S2  = S2  + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(VELY_VAR,  ii,jj,kk)
                 S3  = S3  + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(VELZ_VAR,  ii,jj,kk)
                 E   = E   + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(ENER_VAR,  ii,jj,kk)
                 Nel = Nel + ReconstructionMatrix(iNodeX,n)*solnData(DENS_VAR,ii,jj,kk)*solnData(YE_MSCALAR,ii,jj,kk)
              end do

              uCF(iNodeX,ic,jc,kc,iCF_D ) &
              = D   * conv_dens
              uCF(iNodeX,ic,jc,kc,iCF_S1) &
              = S1  * conv_mom * uGF(iNodeX,ic,jc,kc,iGF_h_1)
              uCF(iNodeX,ic,jc,kc,iCF_S2) &
              = S2  * conv_mom * uGF(iNodeX,ic,jc,kc,iGF_h_2)
              uCF(iNodeX,ic,jc,kc,iCF_S3) &
              = S3  * conv_mom * uGF(iNodeX,ic,jc,kc,iGF_h_3)
              uCF(iNodeX,ic,jc,kc,iCF_E ) &
              = E   * conv_enr
              uCF(iNodeX,ic,jc,kc,iCF_Ne) &
              = Nel * conv_ne

           end do
        end do
     end do
  end do

  ! Copy radiation vars from Flash arrays into Thornado arrays
#if defined(THORNADO_OMP_OL)
  !$OMP TARGET TEAMS DISTRIBUTE COLLAPSE(7)
#elif defined(THORNADO_OACC)
  !$ACC PARALLEL LOOP GANG COLLAPSE(7) &
  !$ACC PRIVATE( i, j, k, ii, jj, kk, ioff, ivar, iNodeE, iNodeX ) &
  !$ACC PRESENT( solnData, uCR, u_lo, u_hi, lo, hi )
#elif defined(THORNADO_OMP)
  !$OMP PARALLEL DO COLLAPSE(7)
#endif
  do kc = u_lo(KAXIS), u_hi(KAXIS)
     do jc = u_lo(JAXIS), u_hi(JAXIS)
        do ic = u_lo(IAXIS), u_hi(IAXIS)
           do iS = 1, THORNADO_NSPECIES ; do iCR = 1, THORNADO_NMOMENTS
              do iE = 1-THORNADO_SWE, THORNADO_NE+THORNADO_SWE
              do iNode = 1, THORNADO_RAD_NDOF

                 i = lo(IAXIS) + THORNADO_NNODESX*(ic-1)
                 j = lo(JAXIS) + THORNADO_NNODESX*(jc-1)
                 k = lo(KAXIS) + THORNADO_NNODESX*(kc-1)

                 ioff = THORNADO_BEGIN &
                    + (iS -1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)*THORNADO_NMOMENTS) &
                    + (iCR-1)*(THORNADO_NNODESE*(THORNADO_NE+2*THORNADO_SWE)) &
                    + (iE -1 + THORNADO_SWE)*(THORNADO_NNODESE)

                 iNodeE = mod((iNode -1)                 ,THORNADO_NNODESE   ) + 1
                 iNodeX = mod((iNode -1)/THORNADO_NNODESE,THORNADO_FLUID_NDOF) + 1

                 ii     = mod((iNodeX-1)                    ,THORNADO_NNODESX) + i
                 jj     = mod((iNodeX-1)/THORNADO_NNODESX   ,THORNADO_NNODESX) + j
                 kk     = mod((iNodeX-1)/THORNADO_NNODESX**2,THORNADO_NNODESX) + k

                 ivar = ioff + iNodeE - 1

                 uCR(iNode,iE,ic,jc,kc,iCR,iS) = solnData(ivar,ii,jj,kk)

              end do
              end do
           end do ; end do
        end do
     end do
  end do

  return
end subroutine rt_tm_reconstruction
