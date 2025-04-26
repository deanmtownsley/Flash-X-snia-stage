!!***if* source/physics/Multiphase/MultiphaseMain/Multiphase_setFluidProps
!! NOTICE
!!  Copyright 2022 UChicago Argonne, LLC and contributors
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
!!
!!
!!***
!!REORDER(4): solnData, face[xyz]Data

#include "Simulation.h"
#include "constants.h"
#include "Multiphase.h"

subroutine Multiphase_setFluidProps(solnData, facexData, faceyData, facezData, del,&
              logc, higc)

   use Multiphase_data
   use Stencils_interface, ONLY: Stencils_lsCenterProps, &
                                 Stencils_lsFaceProps2d, &
                                 Stencils_lsFaceProps3d, &
                                 Stencils_lsCurvature2d, &
                                 Stencils_lsCurvature3d, &
                                 Stencils_lsNormals2d, &
                                 Stencils_lsNormals3d
   use Timers_interface, ONLY: Timers_start, Timers_stop

!---------------------------------------------------------------------------------------------
   implicit none
   real, pointer, dimension(:, :, :, :) :: solnData, facexData, faceyData, facezData
   real,dimension(MDIM),intent(IN) :: del
   integer, dimension(MDIM),intent(IN) :: logc, higc
   
   integer :: ierr, i, j, k
   real minCellDiag
!-------------------------------------------------------------------------------------------------
   

   call Timers_start("Multiphase_setFluidProps")


   minCellDiag = SQRT(del(DIR_X)**2.+del(DIR_Y)**2.+del(DIR_Z)**2.)

   call Stencils_lsCenterProps(solnData(DFUN_VAR, :, :, :), &
                               solnData(SMHV_VAR, :, :, :), &
                               1., &
                               logc( IAXIS), higc( IAXIS), &
                               logc( JAXIS), higc( JAXIS), &
                               logc( KAXIS), higc( KAXIS), &
                               iSmear=mph_iPropSmear*minCellDiag)

   call Stencils_lsCenterProps(solnData(DFUN_VAR, :, :, :), &
                               solnData(PFUN_VAR, :, :, :), &
                               1., &
                               logc( IAXIS), higc( IAXIS), &
                               logc( JAXIS), higc( JAXIS), &
                               logc( KAXIS), higc( KAXIS))

   solnData(mph_iMuCVar, :, :, :) = solnData(SMHV_VAR, :, :, :)*mph_muGas + &
                                    (1 - solnData(SMHV_VAR, :, :, :))*solnData(mph_iMucVar, :, :, :)

   solnData(mph_iRhoCVar, :, :, :) = solnData(SMHV_VAR, :, :, :)*(1./mph_rhoGas) + &
                                     (1 - solnData(SMHV_VAR, :, :, :))*solnData(mph_iRhoCVar, :, :, :)

#if NDIM < MDIM
   call Stencils_lsFaceProps2d(solnData(DFUN_VAR, :, :, :), &
                               facexData(mph_iRhoFVar, :, :, :), &
                               faceyData(mph_iRhoFVar, :, :, :), &
                               1./mph_rhoGas, &
                               logc( IAXIS), higc( IAXIS), &
                               logc( JAXIS), higc( JAXIS))!, &
                               !iSmear=mph_iPropSmear*minCellDiag)

   call Stencils_lsNormals2d(solnData(RHOC_VAR, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), &
                             logc(IAXIS), higc(IAXIS), &
                             logc(JAXIS), higc(JAXIS))

   call Stencils_lsCurvature2d(solnData(CURV_VAR, :, :, :), &
                               solnData(RHOC_VAR, :, :, :), &
                               solnData(NRMX_VAR, :, :, :), &
                               solnData(NRMY_VAR, :, :, :), &
                               del(DIR_X), del(DIR_Y), &
                               logc(IAXIS), higc(IAXIS), &
                               logc(JAXIS), higc(JAXIS))

#else
   call Stencils_lsFaceProps3d(solnData(DFUN_VAR, :, :, :), &
                               facexData(mph_iRhoFVar, :, :, :), &
                               faceyData(mph_iRhoFVar, :, :, :), &
                               facezData(mph_iRhoFVar, :, :, :), &
                               1./mph_rhoGas, &
                               logc( IAXIS), higc( IAXIS), &
                               logc( JAXIS), higc( JAXIS), &
                               logc( KAXIS), higc( KAXIS))!, &
                               !iSmear=mph_iPropSmear*minCellDiag)

   call Stencils_lsNormals3d(solnData(RHOC_VAR, :, :, :), &
                             solnData(NRMX_VAR, :, :, :), &
                             solnData(NRMY_VAR, :, :, :), &
                             solnData(NRMZ_VAR, :, :, :), &
                             del(DIR_X), del(DIR_Y), del(DIR_Z), &
                             logc(IAXIS), higc(IAXIS), &
                             logc(JAXIS), higc(JAXIS), &
                             logc(KAXIS), higc(KAXIS))

   call Stencils_lsCurvature3d(solnData(CURV_VAR, :, :, :), &
                               solnData(RHOC_VAR, :, :, :), &
                               solnData(NRMX_VAR, :, :, :), &
                               solnData(NRMY_VAR, :, :, :), &
                               solnData(NRMZ_VAR, :, :, :), &
                               del(DIR_X), del(DIR_Y), del(DIR_Z), &
                               logc(IAXIS), higc(IAXIS), &
                               logc(JAXIS), higc(JAXIS), &
                               logc(KAXIS), higc(KAXIS))
#endif

   call Timers_stop("Multiphase_setFluidProps")

   return
end subroutine Multiphase_setFluidProps
