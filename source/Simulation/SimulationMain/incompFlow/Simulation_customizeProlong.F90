!!****if* source/Simulation/SimulationMain/incompFlow/Simulation_customizeProlong
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
!! NAME
!!
!!  Simulation_customizeProlong
!!
!! SYNOPSIS
!!
!!  call Simulation_customizeProlong(integer(IN)  :: beforeOrAfter)
!!
!! DESCRIPTION
!!
!!  The Simulation_customizeProlong interface provides a way to
!!  customize the prolongation of Grid data that normally happens
!!  after an AMR Grid has changed - in particular, the interpolation
!!  of data into blocks that were newly created by refining existing
!!  blocks.
!!
!!  The interface is called twice for each time that the global
!!  prolongation operation is applied: once just before prolongation
!!  gets applied, and then again after prolongation is done.
!!  The single arguments is used to distinguish the calls.
!!
!!  After the refinement, apply the user-defined routine for a 
!!  special prolongation routine from the coarse parents to 
!!  the newly created child blocks.
!!
!!  Currently there are two divergence-free prolongation algorithms available:
!!  First is to use the "injection method" in which the face-centered 
!!  divergence-free magnetic fields of the parent blocks are directly 
!!  injected onto the newly born children blocks without any spatial 
!!  interpolation, therefore, automatically maintain divergence-free
!!  properties on the refined blocks. This method works fine in many cases,
!!  especially for smoothly varying fields. It is available for both 
!!  2D and 3D calculation.
!!  Second is the "Balsara's method" that uses bilinear polynomials that
!!  prolong the divergence-free facevars of the parent block data
!!  to the newly born favevars of the children block data, keeping the
!!  divergence-free properties.
!!  See Paramesh's users guide for more detail.
!!
!! ARGUMENTS
!!
!!   beforeOrAfter : BEFORE when called before prolongation,
!!                   AFTER when called after prolongation.
!!
!! NOTES
!!
!!   The constants BEFORE and AFTER are defined in constants.h.
!!
!!   As of FLASH 3.1.1, a non-stub implementation for use by
!!   MHD simulations is provided in the directory tree location
!!   source/Simulation/SimulationMain/magnetoHD . All simulations
!!   placed under the magnetoHD directory will therefore this
!!   implementation by default.  This is usually desired for
!!   simulations using MHD.  For this reason, simulations using
!!   MHD should be placed under the magnetoHD directory.
!!
!!***

subroutine Simulation_customizeProlong(beforeOrAfter)

#include "Simulation.h"

#ifdef FLASH_GRID_PARAMESH
#if NFACE_VARS > 0
  use physicaldata, ONLY : interp_mask_facex,interp_mask_facey,interp_mask_facez,i_divf_fc_vars
  use paramesh_dimensions, ONLY : nfacevar,nfield_divf
  use IncompNS_data, ONLY : ins_prol_method
  use gr_pmDivpres_mod, ONLY : prol_fc_divpres_init
#endif
#endif

  use Driver_interface,ONLY : Driver_abort

  implicit none
#include "constants.h"
#include "IncompNS.h"

  integer, intent(IN) :: beforeOrAfter ! BEFORE for before, AFTER for after

#ifdef FLASH_GRID_PARAMESH
#if NFACE_VARS > 0
  integer, parameter :: totFaces = 3
  integer, parameter  :: tot_divpres_prol = 2 !(VELC and HVN1)
  integer, dimension(NFACE_VARS),save :: interp_mask_facex_old,&
       interp_mask_facey_old,&
       interp_mask_facez_old
  integer, dimension(totFaces,tot_divpres_prol) :: ins_divpres_fc_vars
  integer :: i
#endif
#endif

  if (beforeOrAfter == BEFORE) then

#ifdef FLASH_GRID_PARAMESH
#if NFACE_VARS > 0
     interp_mask_facex_old(:)=interp_mask_facex(:)
     interp_mask_facey_old(:)=interp_mask_facey(:)
     interp_mask_facez_old(:)=interp_mask_facez(:)

     if (ins_prol_method == INJECTION_PROL) then
        !! zeroth order
        interp_mask_facex(VELC_FACE_VAR)=0
        interp_mask_facey(VELC_FACE_VAR)=0
        interp_mask_facez(VELC_FACE_VAR)=0

        interp_mask_facex(HVN1_FACE_VAR)=0
        interp_mask_facey(HVN1_FACE_VAR)=0
        interp_mask_facez(HVN1_FACE_VAR)=0


     elseif (ins_prol_method == LINEAR_PROL) then
        !! Linear
        interp_mask_facex(VELC_FACE_VAR)=1
        interp_mask_facey(VELC_FACE_VAR)=1
        interp_mask_facez(VELC_FACE_VAR)=1

        interp_mask_facex(HVN1_FACE_VAR)=1
        interp_mask_facey(HVN1_FACE_VAR)=1
        interp_mask_facez(HVN1_FACE_VAR)=1


     elseif (ins_prol_method == QUADRATIC_PROL) then
        !! Quadratic
        interp_mask_facex(VELC_FACE_VAR)=2
        interp_mask_facey(VELC_FACE_VAR)=2
        interp_mask_facez(VELC_FACE_VAR)=2

        interp_mask_facex(HVN1_FACE_VAR)=2
        interp_mask_facey(HVN1_FACE_VAR)=2
        interp_mask_facez(HVN1_FACE_VAR)=2


     elseif (ins_prol_method == DIVPRES_PROL) then

        ! Divergence preserving prolongation: 
        ins_divpres_fc_vars(1:totFaces,1) = VELC_FACE_VAR
        ins_divpres_fc_vars(1:totFaces,2) = HVN1_FACE_VAR

        !Divergence cleaning in block boundaries neighboring old blocks:
        i_divf_fc_vars(1:totFaces,1) = VELC_FACE_VAR
        i_divf_fc_vars(1:totFaces,2) = HVN1_FACE_VAR

        call prol_fc_divpres_init(tot_divpres_prol,totFaces, ins_divpres_fc_vars)

     else
        call Driver_abort&
             ("[Grid_updateRefinement]: unknown prolongation algorithm for face-centered magnetic fields!")
     endif
#endif  
#endif

  else if (beforeOrAfter == AFTER) then

#ifdef FLASH_GRID_PARAMESH
#if NFACE_VARS > 0
     interp_mask_facex(:)=interp_mask_facex_old(:)
     interp_mask_facey(:)=interp_mask_facey_old(:)
     interp_mask_facez(:)=interp_mask_facez_old(:)
#endif
#endif


  else
        
     call Driver_abort("Simulation_customizeProlong: this is meaningless!")
  end if
     

end subroutine Simulation_customizeProlong
