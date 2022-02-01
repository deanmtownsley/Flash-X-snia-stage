!!****f* source/Grid/Grid_smoothVar
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
!!  NAME 
!!
!!  Grid_smoothVar
!!
!!  SYNOPSIS
!!
!! 
!!  call Grid_smoothVar(integer(in) :: ivar, 
!!                           integer(in) :: ivarOut,
!!                           real(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:),
!!                           integer(in),value :: lbIU,lbUJ,lbUK,
!!                  OPTIONAL,integer(in) :: smoothMethod,
!!                  OPTIONAL,integer(IN) :: gcLayers,
!!                  OPTIONAL,integer(IN) :: blockID ,
!!                  OPTIONAL,logical(IN) :: useMinSmoothVarVal,
!!                  OPTIONAL,real(IN)    :: minSmoothVarVal,
!!                  OPTIONAL,logical(IN) :: useMaxSmoothVarVal,
!!                  OPTIONAL,real(IN)    :: maxSmoothVarVal,
!!                  OPTIONAL,real(IN)    :: smoothCoeff )
!!
!!  DESCRIPTION 
!!      Smooths one variable in the data passsed in as array solnData.
!!
!! ARGUMENTS
!!
!!   ivar  : index into solution vector, indicating a variable to be smoothed.
!!   ivarOut: index into solution vector, indicating where the smoothed
!!            variable is to be stored.
!!            May be the same as ivar.
!!   solnData : The block data on which to operate
!!   lbUI,lbUJ,lbUK: Lower bounds of assumed-shape array solnData
!!   smoothMethod  : Method for smoothing.
!!   blockID  : Identifies the block on to whic hdata belongs; optional,
!!              for debugging.
!!   gcLayers : In how many layers of guard cells, in addition to the
!!              interior cells, the smoothed output is desired.
!!              Note that one more layer is required in the variable given
!!              by ivar for input to the smoothing.
!!
!! SIDE EFFECTS
!!
!!  Modifies a variable (as requested by argument ivarOut) in the solution
!!  storage UNK.
!!
!!***
subroutine Grid_smoothVar(ivar, ivarOut, &
     solnData, lbUI,lbUJ,lbUK, &
     blklim, smoothMethod, gcLayers, blockID,&
     useMinSmoothVarVal,&
     minSmoothVarVal,&
     useMaxSmoothVarVal,&
     maxSmoothVarVal,&
     smoothCoeff )

  implicit none

#include "constants.h"
#include "FortranLangFeatures.fh"
  
  ! Arguments:
  integer, intent(in) :: ivar
  integer, intent(in) :: ivarOut
  integer, VALUE_INTENT(IN) :: lbUI,lbUJ,lbUK
  real,    intent(INOUT) :: solnData(:,lbUI:,lbUJ:,lbUK:)
  integer, intent(in)    :: blklim(2,MDIM)
  integer, intent(in),OPTIONAL :: smoothMethod
  integer, intent(IN),OPTIONAL :: gcLayers
  integer, intent(IN),OPTIONAL :: blockID
  logical, intent(IN),OPTIONAL :: useMinSmoothVarVal,useMaxSmoothVarVal
  real   , intent(IN),OPTIONAL :: minSmoothVarVal,maxSmoothVarVal
  real   , intent(IN),OPTIONAL :: smoothCoeff

end subroutine Grid_smoothVar
