!!
!! Dean M. Townsley 2009-2026
!!
!! Replacement of default refinement marking from
!!  source/Grid/Grid_markRefineDerefine.F90
!! see that file for API spec and notes.
!!
!! The main difference between the default refinement marking and this
!! implementation is the usage of several different maximum refinement levels
!! which depend on the local physical variables (density, energy
!! generation rate, abundances).  This requires that the parent-child
!! reconciliation (i.e. a child should not derefine if its parent is marked
!! for refinement) be performed after all refinement checks have been
!! performed.  This is somewhat complex because parents may be on different
!! processors than their children, so refinement information must be
!! communicated.
!!
!! The procedure is thus:
!!   1. Use standard routines to mark based on gradients
!!   2. Impose additional limits on refinement level
!!      2.1. impose specified max refinement in fluff
!!         where fluff is determined by flff > FluffThresh (typically 0.99)
!!         regions with minflff > threshold (typcially 0.99) are derefined
!!      2.2. impose specified max refinement outside energy-generating regions
!!         regions are determined by an enuc threshold
!!         note that the margin is on the other side of the
!!         threshold compared to the Fluff threshold above
!!         i.e. regions with enuc > threshold can fully refine
!!      2.3. standard refinement min and max
!!      2.4. ejecta phase resolution limits
!!   3. Do parent-child consistency
!!
!!***

subroutine Grid_markRefineDerefine()

  use Grid_data, ONLY : gr_refine_cutoff, gr_derefine_cutoff,&
                        gr_refine_filter,&
                        gr_numRefineVars,gr_refine_var, &
                        gr_meshComm, gr_meshMe
  use tree, ONLY : newchild, refine, derefine, stay, lrefine, &
                   lrefine_max, lrefine_min, &
                   nodetype, parent, child, nchild, lnblocks
  use Grid_interface, ONLY : Grid_fillGuardCells, &
       Grid_getCellCoords, Grid_getTileIterator, Grid_releaseTileIterator
  use gr_interface,   ONLY : gr_markRefineDerefine
  use Simulation_data, ONLY : sim_refFluffThresh, sim_refFluffMargin, &
       sim_refFluffLevel, sim_refNogenEnucThresh, sim_refNogenMargin, &
       sim_refNogenLevel, sim_refBurnedKeyProductIndex, sim_refBurnedProductThresh, &
       sim_refEjectaPhaseStartTime, sim_refEjectaPhaseMaxRes
  use Driver_interface, ONLY: Driver_getSimTime
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile, ONLY : Grid_tile_t

#include "Flashx_mpi_implicitNone.fh"
#include "Simulation.h"
#include "constants.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,iref

  logical :: gcMask(NUNK_VARS)
  
  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  type(Grid_iterator_t) :: itor
  type(Grid_tile_t) :: tile
  integer ,dimension(MAXBLOCKS) :: blkList
  integer ::  i, j, k

  real, dimension(:,:,:,:), pointer :: solnData

  real :: maxenuc, minflff

  ! for max radius of burned material
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  real :: radiussq, maxrsqburned, global_maxrsqburned, maxRadiusBurned

  real :: simTime, cellsize, ejecta_min_cellsize
  real, dimension(MDIM) :: deltas

  real, dimension(MAXBLOCKS) :: err

  ! for reconciling parents and children
  logical, dimension(MAXBLOCKS) :: refine_parent
  integer :: nsend, nrecv, ierr
  integer, dimension(MAXBLOCKS) :: recvreq
  integer, dimension(MAXBLOCKS*nchild) :: sendreq
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS) :: recvstat
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS*nchild) :: sendstat

  call Driver_getSimTime( simTime )

 !----------------------------------------------------
  !  1. Use standard routines to mark based on gradients
  !----------------------------------------------------

  gcMask = .FALSE.
  do i = 1,gr_numRefineVars
     iref = gr_refine_var(i)
     if (iref > 0) gcMask(iref) = .TRUE.
  end do

  ! refinement max/min levels are based on these variables
  gcMask(FLFF_MSCALAR) = .true.
  gcMask(ENUC_VAR) = .true.

  call Grid_fillGuardCells( CENTER, ALLDIR, eosMode=MODE_DENS_EI, &
       doEos=.true., maskSize=NUNK_VARS, mask=gcMask, makeMaskConsistent=.true.)

  ! first find max radius of burned products, used to set refinement in
  ! ejecta phase
  call Grid_getTileIterator( itor, LEAF )
  maxrsqburned = 0.0
  do while( itor%isValid() )

     ! get info about block
     call itor%currentTile( tile )
     call tile%getDataPtr( solnData, CENTER )

     allocate( iCoords( tile%grownLimits(LOW,IAXIS) : tile%grownLimits(HIGH,IAXIS) ) )
     allocate( jCoords( tile%grownLimits(LOW,JAXIS) : tile%grownLimits(HIGH,JAXIS) ) )
     allocate( kCoords( tile%grownLimits(LOW,KAXIS) : tile%grownLimits(HIGH,KAXIS) ) )
     call Grid_getCellCoords( IAXIS, CENTER, tile%level, tile%grownLimits(LOW,:), tile%grownLimits(HIGH,:), iCoords)
     call Grid_getCellCoords( JAXIS, CENTER, tile%level, tile%grownLimits(LOW,:), tile%grownLimits(HIGH,:), jCoords)
     call Grid_getCellCoords( KAXIS, CENTER, tile%level, tile%grownLimits(LOW,:), tile%grownLimits(HIGH,:), kCoords)

     ! loop over all zones to find max r**2 for burned material
     ! burned is determined by key product above specified threshold
     do k = tile%limits(LOW,KAXIS), tile%limits(HIGH,KAXIS)
        do j = tile%limits(LOW,JAXIS), tile%limits(HIGH,JAXIS)
           do i = tile%limits(LOW,IAXIS), tile%limits(HIGH,IAXIS)
              if ( solnData(sim_refBurnedKeyProductIndex,i,j,k) >= sim_refBurnedProductThresh ) then
                 radiussq = iCoords(i)**2
                 if (NDIM >= 2) radiussq = radiussq + jCoords(j)**2
                 if (NDIM == 3) radiussq = radiussq + kCoords(k)**2
                 maxrsqburned = max( maxrsqburned, radiussq )
              endif
           end do
        end do
     end do
     deallocate(iCoords)
     deallocate(jCoords)
     deallocate(kCoords)

     call tile%releaseDataPtr( solnData, CENTER )

     call itor%next()
  end do
  call Grid_releaseTileIterator( itor )

  call MPI_Allreduce( maxrsqburned, global_maxrsqburned, 1, FLASH_REAL, MPI_MAX, gr_meshComm, ierr)

  maxRadiusBurned = sqrt(global_maxrsqburned)

  ! if (mesh_mype==MASTER_PE) print*, 'max radius burned material:', maxRadiusBurned


  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

  ! default is to derefine unless a criteria requires resolution
  ! note only mark active blocks, otherwise nothing works
  call Grid_getTileIterator( itor, ACTIVE_BLKS, tiling=.false. )
  do while (itor%isValid())
     call itor%currentTile( tile )
     derefine(tile%id) = .true.
     call itor%next()
  enddo
  call Grid_releaseTileIterator( itor )
  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     err(:)      = 0.0
     call gr_estimateError(err, iref, ref_filter)
     call gr_markRefineDerefine(err, ref_cut, deref_cut)
  end do

  !----------------------------------------------------
  !  2. Impose additional limits on refinement level
  !----------------------------------------------------
  call Grid_getTileIterator( itor, ACTIVE_BLKS, tiling=.false. )
  do while (itor%isValid())
     ! get info about block
     call itor%currentTile( tile )
     call tile%getDataPtr( solnData, CENTER )

     ! find max enuc and min flff in block + guardcells
     maxenuc = 0.0
     minflff = 1.0
     do k = tile%grownLimits(LOW,KAXIS), tile%grownLimits(HIGH,KAXIS)
        do j = tile%grownLimits(LOW,JAXIS), tile%grownLimits(HIGH,JAXIS)
           do i = tile%grownLimits(LOW,IAXIS), tile%grownLimits(HIGH,IAXIS)
              maxenuc = max(abs(solnData(ENUC_VAR,i,j,k)),maxenuc)
              minflff  = min(solnData(FLFF_MSCALAR,i,j,k),minflff)
           enddo
        enddo
     enddo

     call tile%releaseDataPtr( solnData, CENTER )

     !----------------------------------------------------
     !  2.1. impose specified max refinement in fluff
     !     where fluff is determined by flff > FluffThresh (typically 0.99)
     !  regions with minflff > threshold (typcially 0.99) are derefined
     !----------------------------------------------------

     ! check min flff in block (including guardcells)
     !  i.e. block must contain only fluff in order to be prevented from refining
     ! above threshold force derefinement to sim_refFluffLevel
     if ( minflff > sim_refFluffThresh ) then
        if ( lrefine(tile%id) > sim_refFluffLevel ) then
           refine(tile%id) = .false.
           derefine(tile%id) = .true.
        else if ( lrefine(tile%id) == sim_refFluffLevel ) then
           refine(tile%id) = .false.
        endif
     ! below threshold but within margin prevent refinement past FluffLevel
     ! but don't force derefinement
     else if ( minflff > sim_refFluffThresh*(1.0-sim_refFluffMargin) &
               .and. lrefine(tile%id) >= sim_refFluffLevel ) then
        refine(tile%id) = .false.
     endif
     ! else allow refinement (lrefine_max is enforced below)

     !----------------------------------------------------
     !  2.2. impose specified max refinement outside energy-generating regions
     !     regions are determined by an enuc threshold
     !  note that the margin is on the other side of the
     !  threshold compared to the Fluff threshold above
     !   i.e. regions with enuc > threshold can fully refine
     !----------------------------------------------------

     ! below threshold and margin force derefinement to NogenLevel
     if ( maxenuc < sim_refNogenEnucThresh*(1.0-sim_refNogenMargin) ) then
        if ( lrefine(tile%id) > sim_refNogenLevel ) then
           refine(tile%id) = .false.
           derefine(tile%id) = .true.
        else if ( lrefine(tile%id) == sim_refNogenLevel ) then
           refine(tile%id) = .false.
        endif
     ! below threshold but within margin prevent refinement past NogenLevel
     ! but don't force derefinement
     else if ( maxenuc < sim_refNogenEnucThresh &
               .and. lrefine(tile%id) >= sim_refNogenLevel ) then
        refine(tile%id) = .false.
     endif
     ! else allow refinement (lrefine_max is enforced below)

     !----------------------------------------------------
     !  2.3. standard refinement min and max
     !----------------------------------------------------

     if (lrefine(tile%id) < lrefine_min) then
        refine(tile%id) = .true.
        derefine(tile%id) = .false.
     else if (lrefine(tile%id) == lrefine_min) then
        derefine(tile%id) = .false.
     endif
     if (lrefine(tile%id) > lrefine_max) then
        refine(tile%id) = .false.
        derefine(tile%id) = .true.
     else if (lrefine(tile%id) == lrefine_max) then
        refine(tile%id) = .false.
     endif
     
     !----------------------------------------------------
     !  2.4. ejecta phase resolution limits
     !----------------------------------------------------

     if (simTime > sim_refEjectaPhaseStartTime) then

        ejecta_min_cellsize = maxRadiusBurned / sim_refEjectaPhaseMaxRes
        call tile%deltas( deltas )
        cellsize = deltas(1)

        ! if cell size for this block is less than min, force derefinement
        if ( cellsize < ejecta_min_cellsize ) then
           derefine(tile%id) = .true.
           refine(tile%id) = .false.
        ! only allow refinement if new cell size would be 110% or larger than minimum
        else if ( 0.5*cellsize < 1.1*ejecta_min_cellsize ) then
              refine(tile%id) = .false.
        endif

     endif

     call itor%next()

  enddo

  call Grid_releaseTileIterator( itor )


  !---------------------------------------------------------------
  ! 3. Do parent-child consistency
  !
  ! for children that are marked derefine, check if parent is 
  ! marked refine and if so unmark derefine
  !---------------------------------------------------------------
  ! just using raw paramesh tree structures here
  refine_parent(1:lnblocks) = .false.
  ! open (async) message recieve if parent is off-procssor
  ! otherwise fill directly
  !    message id is child block number on local processor
  nrecv = 0
  do i = 1, lnblocks
     if (parent(1,i) > -1) then
        if (parent(2,i)/=gr_meshMe) then
           nrecv = nrecv+1
           call MPI_IRecv(refine_parent(i), 1, MPI_LOGICAL, parent(2,i), &
                                    i, gr_meshComm, recvreq(nrecv), ierr)
        else
           refine_parent(i) = refine(parent(1,i))
        endif
     endif
  end do
  ! parents send refine flag to each off-processor child
  nsend = 0
  do i = 1, lnblocks
     do j = 1,nchild
        if (child(1,j,i) > -1) then
           if (child(2,j,i) /= gr_meshMe) then
              nsend = nsend + 1
              call MPI_ISend(refine(i), 1, MPI_LOGICAL, child(2,j,i), &
                                  child(1,j,i), gr_meshComm, sendreq(nsend), ierr)
           endif
        endif
     enddo
  enddo
  ! wait to recieve all parent info
  if (nrecv > 0) then
     call MPI_Waitall(nrecv,recvreq,recvstat,ierr)
  endif
  ! now reconcile, deferring to parent marked for refine
  do i = 1, lnblocks
     if (nodetype(i) == LEAF .and. derefine(i) .and. refine_parent(i) ) then
        derefine(i) = .false.
     endif
  enddo
  ! wait until last to ask to have refine() buffer back under our control since
  ! we are not modifying it
  if (nsend > 0) then
     call MPI_Waitall(nsend,sendreq,sendstat,ierr)
  endif


  ! When the flag arrays are passed to Paramesh for processing, only leaf
  ! blocks should be marked. (noted in default version of this subroutine)
  where (nodetype(:) .NE. LEAF)
     refine(:)   = .false.
     derefine(:) = .false.
  end where


  return
end subroutine Grid_markRefineDerefine
