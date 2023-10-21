!!****f* source/Grid/Grid_fillGuardCells
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
!!  Grid_fillGuardCells
!!
!! SYNOPSIS
!!
!!  call Grid_fillGuardCells(integer(IN) :: gridDataStruct,
!!                           integer(IN) :: idir,
!!                  optional,integer(IN) :: minLayers,
!!                  optional,integer(IN) :: eosMode,
!!                  optional,logical(IN) :: doEos,
!!                  optional,integer(IN) :: maskSize,
!!                  optional,logical(IN),dimension(:) :: mask(maskSize),
!!                  optional,logical(IN) :: makeMaskConsistent,
!!                  optional,logical(IN) :: doLogMask,
!!                  optional,integer(IN) :: selectBlockType,
!!                  optional,logical(IN) :: unitReadsMeshDataOnly)
!!
!!
!! DESCRIPTION
!!
!!  This routine fills the guardcells of the specified data structure of 
!!  each block present in the mesh. If the adjacent blocks are at the same
!!  level of refinement, then the values are copied from the appropriate
!!  internal cells of the neighboring block. If the blocks are at different
!!  levels of refinement, then in addition to copying, there will be 
!!  restriction or prolongation. 
!!
!!  The optional arguments related to the Eos calls are provided because if
!!  a solver relies on a call to Eos to establish thermodynamic equilibrium
!!  in various variables, then either the solver can call Eos explicitly itself,
!!  or more conveniently, it can instruct the GC-fill routine to do it. This
!!  feature is especially useful when mesh quantities are used to derive,
!!  for example, tracer particle attributes. If the user opts to enable doEos,
!!  they will not have to determine whether or where an Eos call is needed, the
!!  code will do it for them. Note that Eos calls via Grid_fillGuardCells will
!!  only be applied to guard cells, and potentially only to guard cell regions
!!  where the data in those guard cells may be the result of either
!!  interpolation or cell averaging because a block boundary coincides with or
!!  touches a fine/coarse boundary.
!!
!!  For performance reasons, users may choose to use masking in filling the
!!  guardcells, which is a feature supported only with Paramesh4 and Amrex. They can
!!  do so by using the optional arguments mask, maskSize, and makeMaskConsistent.
!!  However, users should exercise extreme caution in using masking. Please see
!!  the NOTES section for potential side effects. If masking is present, and
!!  both makeMaskConsistent and doEos are true, the Eos calculation is applied only
!!  if at least one of the variables selected in the mask is potentially an output
!!  of the Eos calculation. A local routine gr_makeMaskConsistent or
!!  gr_makeMaskConsistent_gen handles this processing.
!!
!!  The argument "gridDataStruct" can take one of several valid values to determine
!!  a specific grid data structure (or combination of data structures) on which to
!!  apply the guardcell fill operation. The currently available options are listed
!!  with the arguments. Most calls in a typical Flash-X application will use CENTER
!!  as the option, since it is most common to use cell-centered grid data and fill
!!  guardcells for all such variables. When using AMR with Paramesh4, the
!!  subroutine Grid_fillGuardCells can act on an additional cell-centered data
!!  structure provided by Paramesh (by default with space for only a single
!!  variable) represented by the option "WORK". More specialized applications
!!  may want to use other options.
!!  The user can also choose to fill guard cells either in a single direction,
!!  or all of them. For the Flash-X solvers supplied as of 2023,
!!  guard cell are filled in all directions.
!!
!!  
!!
!!
!! ARGUMENTS 
!!  
!!
!!  gridDataStruct - integer constant, defined in "constants.h", indicating for
!!                   which grid data structures the guard cells need to be filled.
!!
!!                   Paramesh has 5 data structures for grid variables, the first
!!                   four store the physical variables defined on the mesh. The
!!                   fifth one includes by default a single variable.
!!
!!                   unk                cell centered,
!!                   facex,facey,facez  face centered along i,j,k 
!!                                      direction, respectively,
!!                   work               cell centered, single variable.
!!                   
!!                   valid values of gridDataStruct are  
!!                   CENTER             unk only
!!                   WORK               work   (only supported with Paramesh4 Grid)
!!                   FACEX              facex  (only supported with Paramesh4 and UG)
!!                   FACEY              facey  (only supported with Paramesh4 and UG)
!!                   FACEZ              facez  (only supported with Paramesh4 and UG)
!!                   FACES              facex,facey, and facez (up to NDIM)
!!                   CENTER_FACES       unk and facex,facey,facez (up to NDIM)
!!
!!  idir - direction of guardcell fill.  User can specify ALLDIR for all (x,y,z)
!!         directions, or if, for example, the algorithm only does one directional
!!         sweep at a time, then time can be saved by filling only the guardcell
!!         direction that is needed.  A user would pass in the constants IAXIS,
!!         JAXIS, or KAXIS defined in constants.h to fill guardcells in only one
!!         direction.
!!         All layers of guardcells in the given direction(s) are filled if one
!!         of IAXIS, JAXIS, or KAXIS is specified, or if ALLDIR is specified
!!         and minLayers (see below) is not present.
!!
!!  minLayers - minimum number of guardcell layers requested for all directions.
!!              The caller requests at least this many layers of
!!              guardcells to be filled.  If idir is given as one of IAXIS,
!!              JAXIS, or KAXIS, this applies to any directions NOT selected by
!!              idir. On the other hand, if idir is given as ALLDIR, then
!!              minLayers appliers to all directions equally.
!!              If not specified, the default is 0 if idir is given as IAXIS,
!!              JAXIS, or KAXIS, meaning no guardcells need to be filled in for
!!              the directions perpendicular to what idir selects; the default
!!              is NGUARD, the full number of available guardcells, if idir
!!              is ALLDIR.
!!              Note that the caller can specify, using minLayers, how many
!!              layers it needs filled, but the implementation may do
!!              more and actually fill all or some additional layers.
!!              The caller must therefore not rely on some guardcells
!!              remaining unchanged.
!!              
!!              This argument is meaningful with AMR Grid implementations only,
!!              and for the most part ignored by the Amrex Grid implementation
!!              (but don't count on it).
!!
!!   eosMode -  the EOS mode being used by the solver that is calling the 
!!              routine. Some of the valid values are :
!!              MODE_DENS_EI     density and energy input, pressure and temp output
!!              MODE_DENS_PRES   density/pressure input, temperature/energy output
!!              MODE_DENS_TEMP   density/temperature input, pressure/energy output
!!              Other values may be supported by specific Eos implementations.
!!              When this optional argument is not present in a call, the behavior
!!              is the same as when the default Eos mode that is being used by the
!!              Grid unit is specified; the latter will usually be determined by
!!              a runtime parameter such as eosModeInit or eosMode.
!!
!!  doEos    - a logical variable indicating if the calling routine wants the
!!             GC-fill routine to also make sure that Eos is applied to achieve
!!             thermodynamically consistent values of all variables.
!!
!!
!!             The mask-related arguments are useful only with PARAMESH4. They have
!!             no meaning in UG.  The Amrex Grid currently provides a partial
!!             !DEV: implementation with known failures.
!!
!!  maskSize - size of the part of the mask array to be used.
!!             Should normally be the same as the actual size of the mask array.
!! 
!!  mask -  It is a one-dimensional logical array
!!          with indices corresponding to variables in the grid data
!!          structures. If a variable should have its guardcells filled,
!!          the corresponding element in "mask" is true, otherwise it is
!!          false.
!!          The mask is always ignored if the runtime parameter
!!          enableMaskedGCFill is set .FALSE.
!!  
!! makeMaskConsistent - If true then when a mask is applied, the implementation of
!!          Grid_fillGuardCells makes sure that in addition to all the variables
!!          that are explicitly selected by the mask, all other variables (of
!!          the same grid data structure selected by gridDataStruct) on which
!!          those variables depend for proper GC-fill operation are also
!!          included in the set of variables that are communicated and updated
!!          in guard cells. (See NOTES section for an example regarding the
!!          dependence of conserved variables on the density variable.)
!!          the the doEos argument is also present and true, then the
!!          implementation of Grid_fillGuardCells also determines whether there
!!          is a need to apply the Eos, and the Eos calculation may be skipped
!!          if not.
!!
!! doLogMask - If present and true when a mask is applied, stamp some information
!!          about masking (and whether Eos is or would be called) to the Logfile.
!!          Ignored in implementations that do not support variable masking.
!!
!!  selectBlockType - selects the blocks whose guard cells must be filled
!!            by block type.
!!
!!              This argument is ignored in UG Implementations.
!!
!!              For PARAMESH Grid implementations, recognized values are :
!!
!!              ALL_BLKS    all local blocks on a processor.
!!                          This is not valid in all PARAMESH versions.
!!
!!              ACTIVE_BLKS all currently active blocks, in the PARAMESH
!!                          context this means parent and leaf blocks,
!!                          that is, blocks whose node type is 1 or 2,
!!                          but not blocks of node type 3 (ANCESTOR).
!!              
!!              LEAF        only LEAF blocks, i.e., blocks whose PARAMESH
!!                          node type is 1
!!            
!!              These constants are defined in constants.h
!!
!!              When selectBlockType is not present, the default behavior
!!              with a Paramesh4 Grid is essentially as described for
!!              ACTIVE_BLKS. However, note that there may be some blocks
!!              of PARENT type that have some child blocks that are NOT
!!              of LEAF type (but are themselves of PARENT type). Such
!!              PARENT blocks may not get fully updated in regions that
!!              are not covered by LEAF blocks.
!!
!!       Note that if advance_all_levels is set (in a PARAMESH version
!!       that implements this global flag), guard cells in all levels of
!!       blocks are filled anyway.
!!
!!       Note that the caller can specify, using selectBlockType, which
!!       blocks it needs filled, but the implementation may do
!!       more and actually fill guard cells in additional blocks.
!!       The caller must therefore not rely on guardcells
!!       in other blocks remaining unchanged.
!!
!! unitReadsMeshDataOnly - specifies that the unit calling Grid_fillGuardCells
!!                         will not update any internal grid data after the call;
!!                         or, if it does, it will notify the Grid unit of it
!!                         by calling Grid_notifySolnDataUpdate later.
!!                         This protocol allows the NEXT call of
!!                         Grid_fillGuardCells to skip the actual communication
!!                         when it detects that the guard cells already contain
!!                         up to date data.
!!
!! EXAMPLE 
!!
!!   #include "Simulation.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( CENTER, IAXIS)
!!
!!     This call will fill all guardcells for all cell-centered 
!!     variables in the x direction.
!!     
!! EXAMPLE 2
!!
!!   #include "Simulation.h"
!!   #include "constants.h"
!!
!!      call Grid_fillGuardCells( WORK, ALLDIR)
!!     
!!     This call fills guardcells along all directions. The operation is applied
!!     to the WORK data structure available in Paramesh only.
!!
!! SIDE EFFECTS
!!
!!  After this function returns, all parents of leaf blocks will have current and 
!!  valid solution data (at least for the variables determined by the gridDataStruct
!!  and mask dummy arguments). This is because amr_guardcell calls amr_restrict
!!  internally.  If selectBlockType is used, even more parents of blocks of interest
!!  may get updated by restriction.
!!
!! NOTES
!!
!!  In the default mode, this routine fills guardcells of all the variables
!!  in the specified data structure. However, if masking is used
!!  it is possible to fill guardcells of only 
!!  selected variables. The users must exercise a great deal of caution in 
!!  using masking in the guardcell filling, since masking out some variables may
!!  have unintended consequences.  For example, if any conserved variable
!!  is set to true in the mask, density must be true, otherwise the interpolated
!!  values of the the conserved variable at fine-coarse boundaries will be 
!!  wrong. It is highly recommended that the argument makeMaskConsistent be set to 
!!  true when using masking, since that provides error checking for
!!  the variable sets included in the released solvers. However,  
!!  these checks may not be enough for new solvers introduced by the users.
!!
!!  This function, or one of the lower-level functions invoked by it, MUST
!!  be called (with updated solution data) before child blocks may be removed
!!  in the course of derefinement! See side effects, above, for the reason.
!!
!! 
!!***



subroutine Grid_fillGuardCells( gridDataStruct, idir,&
     minLayers, eosMode, doEos, maskSize, mask, makeMaskConsistent,&
     doLogMask,selectBlockType,unitReadsMeshDataOnly)
  
  implicit none
  integer, intent(in) :: gridDataStruct
  integer, intent(in) :: idir
  integer, optional,intent(in) :: minLayers
  integer, optional, intent(in) :: eosMode
  logical, optional,intent(IN) :: doEos
  integer, optional, intent(IN) :: maskSize
  logical, optional, dimension(:),intent(IN) :: mask
  logical, optional,intent(IN) :: makeMaskConsistent
  logical, optional,intent(IN) :: doLogMask
  integer, optional,intent(IN) :: selectBlockType
  logical, optional, intent(IN) :: unitReadsMeshDataOnly
end subroutine Grid_fillGuardCells
