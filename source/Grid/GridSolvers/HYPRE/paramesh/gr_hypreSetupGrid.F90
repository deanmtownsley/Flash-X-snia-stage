!!****if* source/Grid/GridSolvers/HYPRE/paramesh/gr_hypreSetupGrid
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
!!   gr_hypreSetupGrid
!!
!!  SYNOPSIS
!!
!!  call gr_hypreSetupGrid(integer(IN)  :: blockCount,
!!                         integer(IN), dimension(blockCount) :: blockList,
!!                OPTIONAL,integer(IN)  :: nvars)
!!
!!
!!  DESCRIPTION
!!   This routine sets up the HYPRE Grid (multiple calls in paramesh).
!!
!!
!! ARGUMENTS
!!   blockCount     : The number of blocks in the list.
!!   blockList      : The list of blocks on which the solution must be updated.
!!   nvars          : Number of variables, also number of equations, for a
!!                    system. Default is 1.
!!
!! SIDE EFFECTS
!!
!!  HYPRE objects for the handles gr_hypreMatA, gr_hypreVecB, and gr_hypreVecX
!!  are created and initialized. These will be accessible via the SStruct
!!  interface of HYPRE.
!!
!! NOTES
!!
!!   Uses HYPRE library.
!!
!!   Typically called from gr_hypreGridStatus, which detects when
!!   re-setup is necessary in a running simulation.
!!
!!   Arrays gr_hypreNeghLevels and gr_hypreSurrBlkSum are allocated and
!!   defined to stores the refinement levels of neighbor blocks and
!!   surrounding block summary info. Arrays gr_hypreLower and gr_hypreUpper
!!   store lower and upper indices for each direction for inner cells of
!!   each block. (These indices are specific to the HYPRE part to which
!!   a block belongs based on its refinement level.)
!!
!!   Some other global info is stored in module variables of gr_hypreData,
!!   like gr_hypreRefineMAX, gr_hypreRefineMIN, gr_hypreNParts, etc.
!!
!! SEE ALSO
!!
!!  gr_hypreGridStatus
!!
!!***

subroutine gr_hypreSetupGrid (blockCount, blockType, nvars)

  use Grid_interface,   ONLY : Grid_getTileIterator, Grid_releaseTileIterator

  use Grid_tile, ONLY : Grid_tile_t
  use Grid_iterator, ONLY : Grid_iterator_t
 
  use gr_hypreData,     ONLY : gr_hypreGridIsSetUp, gr_hypreLower, gr_hypreUpper, gr_hypreVecX, &
                               gr_hypreVecB, gr_hypreMatA, gr_hypreGraph, gr_hypreStencils, gr_hypreNParts, &
                               gr_hypreNVars, gr_hypreGrid, gr_hypreRefineMIN, gr_hypreRefineMAX, gr_hypreNeghLevels, &
                               gr_hypreSolverType, gr_hypreSurrBlkSum, gr_hypreMatrixIsSetup

  use Grid_data,        ONLY : gr_meshComm, gr_meshMe

  use Driver_interface, ONLY : Driver_abort

  use Timers_interface, ONLY : Timers_start, Timers_stop

  use gr_hypreLocalInterface, ONLY: gr_hypreAddGraph, gr_hypreLevel
!  use gr_interfaceTypeDecl
  use gr_interface,     ONLY : gr_findAllNeghID, gr_getBlkHandle
  use tree,             ONLY:  surr_blks, lrefine_max


  implicit none

#include "Simulation.h"
#include "constants.h"
#include "HYPREf.h"
#include "Flashx_mpi.h"

  integer,                      intent(IN) :: blockCount
  integer,                      intent(IN) :: blockType
  integer,OPTIONAL,             intent(IN) :: nvars

  integer :: offsets(2*NDIM+1,NDIM)
  integer :: datasize(MDIM)
  integer :: lb, blockID
  integer :: ilower(MDIM), iupper(MDIM), to_index(MDIM), from_index(MDIM)
  integer :: mylower(MDIM)
  integer :: stride (MDIM)
  integer :: part
  integer :: object_type
  integer :: ent, var, i,j,k, eachNegh, ii, jj
  integer :: nentries
  integer :: ierr, numNegh
  integer, dimension(2,MDIM) :: blkLimitsGC, blkLimits
  integer :: lrefine(blockCount), gmax, gmin
 
  type(Grid_tile_t) :: tileDesc
  type(Grid_iterator_t) :: itor
 
  integer, dimension(BLKNO:PROCNO) :: neghBlkProc
  integer ::neghBlk

  integer :: allCenters, iface
  integer :: graph2part, level

  !!     This comes from 'sstruct_mv/HYPRE_sstruct_mv.h'
  integer ::   HYPRE_SSTRUCT_VARIABLE_CELL = 0 !! CELL CENTERED

  integer,dimension(BLKNO:TYPENO) :: negh_prop

  !! Do nothing, if grid already setup, just return.
  if (gr_hypreGridIsSetUp) return

  call Timers_start("gr_hypreSetupGrid")

  gr_hypreGridIsSetUp = .TRUE.
  gr_hypreMatrixIsSetup = .FALSE.

  if (gr_hypreSolverType == HYPRE_SPLIT) then
     object_type = HYPRE_SSTRUCT
  else
     object_type = HYPRE_PARCSR
  end if

  !! Lower and upper Corner ID of each local leaf block
  allocate (gr_hypreLower(blockCount, MDIM))
  allocate (gr_hypreUpper(blockCount, MDIM))

  !! Stores the refinement level of neighbor along each direction.
  allocate (gr_hypreNeghLevels(blockCount, 1:RIGHT_EDGE, 1:RIGHT_EDGE, 1:RIGHT_EDGE))

  allocate (gr_hypreSurrBlkSum(blockCount))


  gr_hypreRefineMAX = 0
  gr_hypreRefineMIN = 9999

  !call gr_hypreLevel()

  datasize = 0
  !!-----------------------------------------------------------------------
  !!     1.  Store lower corner ID, upper corner ID of local leaf blocks
  !!     2.  Compute max/min refinement levels on all leaf blocks.
  !!     3.  Store stride of all local leaf blocks.
  !!-----------------------------------------------------------------------
  lb = 1
  call Grid_getTileIterator(itor, nodetype=blockType)
  do while(itor%isValid())

     call itor%currentTile(tileDesc)

     blkLimits   = tileDesc%limits
     blkLimitsGC = tileDesc%blkLimitsGC
     ilower      = tileDesc%cid
     stride      = tileDesc%stride
 
     ilower(1:MDIM) =  ceiling(real(ilower(1:MDIM)) / real(stride(1:MDIM)))
     datasize  (1:MDIM)= max(blkLimits(HIGH,1:MDIM)-blkLimits(LOW,1:MDIM), datasize(1:MDIM))
     iupper(1:MDIM) = ilower(1:MDIM) + datasize (1:MDIM)

     do i=1, NDIM
        gr_hypreLower(lb,i) =  ilower(NDIM-i+1)
        gr_hypreUpper(lb,i) =  iupper(NDIM-i+1)
     end do

     lrefine(lb) = tileDesc%level
     gr_hypreRefineMAX = max(lrefine(lb), gr_hypreRefineMAX)
     gr_hypreRefineMIN = min(lrefine(lb), gr_hypreRefineMIN)

     lb = lb + 1
     call itor%next()

  end do
  call Grid_releaseTileIterator(itor) 


  do i=NDIM+1, MDIM
     gr_hypreLower(:,i) = 1
     gr_hypreUpper(:,i) = 1

  end do

  !! The max/min refinement levels could reside in different blocks.
  call mpi_allreduce (gr_hypreRefineMAX, gmax, 1, FLASH_INTEGER ,MPI_MAX, gr_meshComm, ierr)
  call mpi_allreduce (gr_hypreRefineMIN, gmin, 1, FLASH_INTEGER ,MPI_MIN, gr_meshComm, ierr)

  gr_hypreRefineMAX = gmax
  gr_hypreRefineMIN = gmin

  !!-----------------------------------------------------------------------
  !!     4.  All neighbors at same refinement level in UG.
  !!-----------------------------------------------------------------------



  !!-----------------------------------------------------------------------
  !!     5.  All blocks (box) at a refinement level forms a part.
  !!         NOTE: Total parts in UG would be 1.
  !!     6.  Only one variable to solve.
  !!-----------------------------------------------------------------------
  gr_hypreNParts = gr_hypreRefineMAX - gr_hypreRefineMIN + 1 !!
  if (present(nvars)) then
     gr_hypreNVars  = nvars
  else
     gr_hypreNVars  = 1
  end if

  part   = 0   !! part iterator


  !!-----------------------------------------------------------------------
  !!     7.  Create a HYPRE grid object with computed number of
  !!         parts and given dimension.
  !!-----------------------------------------------------------------------

  call HYPRE_SStructGridCreate(gr_meshComm, NDIM, gr_hypreNParts, gr_hypreGrid, ierr)


  !!-----------------------------------------------------------------------
  !!     9.  Each leaf block in FLASH is a box object (defined by its
  !!         extents) in HYPRE. Add the leaf blocks to their respective parts
  !!         (refinement level in FLASH).
  !!-----------------------------------------------------------------------
  do lb=1, blockCount
     part = lrefine(lb) - gr_hypreRefineMIN
     call HYPRE_SStructGridSetExtents(gr_hypreGrid, part, gr_hypreLower(lb,1:NDIM), gr_hypreUpper(lb,1:NDIM),ierr)
  end do

  !!-----------------------------------------------------------------------
  !!     10. We are solving for cell centered variables,
  !!         let HYPRE know that. Also, set number of variables per part. It
  !!         is possible to have more then one variable per part with
  !!         different variable types.
  !!-----------------------------------------------------------------------

  do part = 0, gr_hypreNParts-1
     call HYPRE_SStructGridSetVariables(gr_hypreGrid, part, gr_hypreNVars, &
          (/(HYPRE_SSTRUCT_VARIABLE_CELL, var=1,gr_hypreNVars)/) ,ierr)
  end do

  !!-----------------------------------------------------------------------
  !!     11. Assemble HYPRE grid (Global call within HYPRE).
  !!-----------------------------------------------------------------------
  call HYPRE_SStructGridAssemble(gr_hypreGrid, ierr)

  !!-----------------------------------------------------------------------
  !!     12. Define the discretization stencil, we use standard central
  !!         differencing for diffusion operator, hence we get 3,5,7 point
  !!         stencil for 1D, 2D and 3D respectively.
  !!----------------------------------------------------------------------
  nentries = NDIM*2 + 1

  !!-----------------------------------------------------------------------
  !!     13. Offsets provide relative positions in Matrix, this saves us
  !!         the effort of computing exact positions on global matrix,
  !!         transformation of offsets to Matrix is done implicitly by HYPRE.
  !!         NOTE:  These offsets are different from the order in the HYPRE
  !!                manual: because we use FORTRAN arrays, the indices
  !!                are reversed.
  !!----------------------------------------------------------------------

#if NDIM == 1
  offsets(1,1) =  0
  offsets(2,1) =  -1
  offsets(3,1) =  +1
#endif

#if NDIM == 2
  offsets(1,1) =  0
  offsets(1,2) =  0

  offsets(2,1) =  0
  offsets(2,2) =  -1

  offsets(3,1) =  0
  offsets(3,2) =  1

  offsets(4,1) =  -1
  offsets(4,2) =  0

  offsets(5,1) =  1
  offsets(5,2) =  0
#endif

#if NDIM == 3
  offsets(1,1) =  0
  offsets(1,2) =  0
  offsets(1,3) =  0

  offsets(2,1) =  0
  offsets(2,2) =  0
  offsets(2,3) =  -1

  offsets(3,1) =  0
  offsets(3,2) =  0
  offsets(3,3) =  1

  offsets(4,1) =  0
  offsets(4,2) =  -1
  offsets(4,3) =  0

  offsets(5,1) =  0
  offsets(5,2) =  1
  offsets(5,3) =  0

  offsets(6,1) =  -1
  offsets(6,2) =  0
  offsets(6,3) =  0

  offsets(7,1) =  1
  offsets(7,2) =  0
  offsets(7,3) =  0
#endif


  !!-----------------------------------------------------------------------
  !!     14. Step done for pure convinience, by assigning specific numerical
  !!         values to offsets, we can easily refer to them later (as when
  !!         needed). In the single-variable case, var is just 0 for our problem.
  !!----------------------------------------------------------------------

  allocate(gr_hypreStencils(0:0))
  call HYPRE_SStructStencilCreate(NDIM, nentries, gr_hypreStencils(0), ierr)

  var    = 0   !! var iterator.
  do ent = 1, nentries
     call HYPRE_SStructStencilSetEntry(gr_hypreStencils(0), ent-1, offsets(ent,1:NDIM), var, ierr)
  enddo

  !!-----------------------------------------------------------------------
  !!     15. Create Graph object, a graph object is a catch all, it is used
  !!         for building up relations across parts, i.e., at fine-coarse
  !!         boundaries, these cannot be handled using normal stenciled objects.
  !!         NOTE: Not needed in UG.
  !!----------------------------------------------------------------------

  call HYPRE_SStructGraphCreate(gr_meshComm, gr_hypreGrid, gr_hypreGraph, ierr)


  !!-----------------------------------------------------------------------
  !!     16. HYPRE provides a variety of storage formats. Specific formats
  !!         support specific solvers. HYPRE_PARCSR supports most of the
  !!         solvers, especially PCG/AMG which we are interested in using.
  !!----------------------------------------------------------------------

  call HYPRE_SStructGraphSetObjectType(gr_hypreGraph,object_type, ierr)


  !!-----------------------------------------------------------------------
  !!     17. Associate a stencil with each part/variable. We use the same
  !!         stencil for all parts.
  !!----------------------------------------------------------------------
  do part = 0, gr_hypreNParts-1
     call HYPRE_SStructGraphSetStencil(gr_hypreGraph, part, var, gr_hypreStencils(0), ierr)
  end do

  !!-----------------------------------------------------------------------
  !!     18. Creating graphs, done if and only if we are using PARAMESH
  !!         and gr_hypreNParts > 1. Which ensure we have a fine-coasrse
  !!         interface.
  !!----------------------------------------------------------------------


  gr_hypreNeghLevels = -1

  if (gr_hypreNParts > 1) then
     !!-----------------------------------------------------------------------
     !!     19. This entire step is devoted to computing the refinement
     !!         levels (i.e HYPRE part number) of our neighbors. The process is
     !!         a bit warped as we need off processor information.
     !!         NOTE: i)  A dependency on PARAMESH is built into the system, as a
     !!                    tradeoff we avoid off processor communication.
     !!               ii) Some block information might be cached by PARAMESH depending
     !!                   on last performed operation.
     !!               iii)We store this information locally.
     !!----------------------------------------------------------------------
     allCenters = 2**NDIM
     lb = 1
     call Grid_getTileIterator(itor, nodetype=blockType)
     do while(itor%isValid())

        call itor%currentTile(tileDesc)
        blockID = tileDesc%id

        call gr_findAllNeghID(blockID, gr_hypreSurrBlkSum(lb))

        do k = LEFT_EDGE , LEFT_EDGE+K3D*(RIGHT_EDGE-1)
           do j = LEFT_EDGE , LEFT_EDGE+K2D*(RIGHT_EDGE-1)
              do i = LEFT_EDGE , RIGHT_EDGE
                 if (allCenters /= i*j*k) then
                    gr_hypreNeghLevels(lb,i,j,k) = lrefine(lb)
                 end if
              end do
           end do
        end do
        lb = lb + 1
        call itor%next()

     end do
     call Grid_releaseTileIterator(itor) 


     !!-----------------------------------------------------------------------
     !!     20. Using the magic of HYPRE graph object each cell in the fine
     !!         coarse boundary or part boundary is associated with its neighbor.
     !!         NOTE: a)Building in a dependency on PARAMESH as Grid_getBlkCornerID
     !!                 can provide information on off-processor blocks.
     !!               b)The graph is done both ways from fine->coarse and
     !!                 coarse to fine.
     !!               c)The process is repeated one side at a time.
     !!----------------------------------------------------------------------
     lb = 1
     call Grid_getTileIterator(itor, nodetype=blockType)
     do while(itor%isValid())

        call itor%currentTile(tileDesc)

        blockID = tileDesc%id
        mylower = tileDesc%cid
        stride  = tileDesc%stride

        part = lrefine(lb) - gr_hypreRefineMIN

        !! BUGFIX 2013-01-07 - KW
        !! For a given cell, the order in which graph entries from that cell are added
        !! here must match the order in which they are given values!
        !! Cf. gr_hypreCreateMatrix.
        do iface = 1, 2*NDIM
           call gr_hypreAddGraph (lb, blockID, part, iface, datasize, mylower,stride)
        end do

        lb = lb + 1
        call itor%next()

     end do
     call Grid_releaseTileIterator(itor) 
  else
     !! PARAMESH in UG mode
     gr_hypreNeghLevels = gr_hypreRefineMIN !! which is same as gr_lrefineMAX
  end if !! nparts > 1

  !!-----------------------------------------------------------------------
  !!     21. Assemble the hYPRE graph object
  !!-----------------------------------------------------------------------
  call HYPRE_SStructGraphAssemble(gr_hypreGraph, ierr)

  !!-----------------------------------------------------------------------
  !!     22. Create empty matrix,vector objects.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixCreate(gr_meshComm, gr_hypreGraph, gr_hypreMatA, ierr)
  call HYPRE_SStructVectorCreate(gr_meshComm, gr_hypreGrid,  gr_hypreVecB, ierr)
  call HYPRE_SStructVectorCreate(gr_meshComm, gr_hypreGrid,  gr_hypreVecX, ierr)


  !!-----------------------------------------------------------------------
  !!     23. Set storage format type.
  !!         As stated before, this format supports PCG/AMG.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixSetObjectTyp(gr_hypreMatA, object_type, ierr)
  call HYPRE_SStructVectorSetObjectTyp(gr_hypreVecB, object_type, ierr)
  call HYPRE_SStructVectorSetObjectTyp(gr_hypreVecX, object_type, ierr)


  !!-----------------------------------------------------------------------
  !!     24. HYPRE allocates memory for objects.
  !!-----------------------------------------------------------------------
  call HYPRE_SStructMatrixInitialize(gr_hypreMatA, ierr)
  call HYPRE_SStructVectorInitialize(gr_hypreVecB, ierr)
  call HYPRE_SStructVectorInitialize(gr_hypreVecX, ierr)

  call Timers_stop("gr_hypreSetupGrid")

end subroutine gr_hypreSetupGrid
