!!
!! Dean M. Townsley 2009
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
                        gr_numRefineVars,gr_refine_var
  use tree, ONLY : newchild, refine, derefine, stay, lrefine, &
                   lrefine_max, lrefine_min, &
                   nodetype, parent, child, nchild
  use Grid_interface, ONLY : Grid_getBlkPtr, Grid_getListOfBlocks,&
       Grid_releaseBlkPtr, Grid_getBlkIndexLimits, Grid_fillGuardCells, &
       Grid_getCellCoords, Grid_getDeltas
  use Simulation_data, ONLY : sim_refFluffThresh, sim_refFluffMargin, &
       sim_refFluffLevel, sim_refNogenEnucThresh, sim_refNogenMargin, &
       sim_refNogenLevel, sim_refBurnedKeyProductIndex, sim_refBurnedProductThresh, &
       sim_refEjectaPhaseStartTime, sim_refEjectaPhaseMaxRes
  use Driver_interface, ONLY: Driver_getMype, Driver_getComm, Driver_getSimTime

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"

  real :: ref_cut,deref_cut,ref_filter
  integer       :: l,i,iref

  logical :: gcMask(NUNK_VARS)
  
  integer, dimension(2,MDIM) :: blkLimits
  integer, dimension(2,MDIM) :: blkLimitsGC
  integer ,dimension(MAXBLOCKS) :: blkList
  integer :: blkCount, iblk, j, k, error, mesh_mype, meshcomm

  real, dimension(:,:,:,:), pointer :: solnData

  real :: maxenuc, minflff

  ! for max radius of burned material
  integer :: isizeGC, jsizeGC, ksizeGC
  real, allocatable, dimension(:) :: iCoords, jCoords, kCoords
  real :: radiussq, maxrsqburned, global_maxrsqburned, maxRadiusBurned

  real :: simTime, cellsize, ejecta_min_cellsize
  real, dimension(MDIM) :: deltas

  ! for reconciling parents and children
  logical, dimension(MAXBLOCKS) :: refine_parent
  integer :: nsend, nrecv, ierr
  integer, dimension(MAXBLOCKS) :: recvreq, sendreq
  integer, dimension(MPI_STATUS_SIZE,MAXBLOCKS) :: recvstat, sendstat

  call Driver_getMype( MESH_COMM, mesh_mype )
  call Driver_getComm( MESH_COMM, meshcomm )
  call Driver_getSimTime( simTime )

  newchild(:) = .FALSE.
  refine(:)   = .FALSE.
  derefine(:) = .FALSE.
  stay(:)     = .FALSE.

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
  call Grid_getListOfBlocks(LEAF, blkList, blkCount)
  maxrsqburned = 0.0
  do iblk = 1, blkCount

     ! get info about block
     call Grid_getBlkPtr( blkList(iblk), solnData, CENTER)
     call Grid_getBlkindexLimits( BlkList(iblk), blkLimits, blkLimitsGC)

     isizeGC = blkLimitsGC(HIGH,IAXIS)
     allocate(iCoords(isizeGC))
     jsizeGC = blkLimitsGC(HIGH,JAXIS)
     allocate(jCoords(jsizeGC))
     ksizeGC = blkLimitsGC(HIGH,KAXIS)
     allocate(kCoords(ksizeGC))
     call Grid_getCellCoords(IAXIS,blkList(iblk),CENTER,.true.,iCoords,isizeGC)
     call Grid_getCellCoords(JAXIS,blkList(iblk),CENTER,.true.,jCoords,jsizeGC)
     call Grid_getCellCoords(KAXIS,blkList(iblk),CENTER,.true.,kCoords,ksizeGC)

     ! loop over all zones to find max r**2 for burned material
     ! burned is determined by key product above specified threshold
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
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
      call Grid_releaseBlkPtr(blkList(iblk),solnData)

  end do

  call MPI_Allreduce( maxrsqburned, global_maxrsqburned, 1, FLASH_REAL, MPI_MAX, meshcomm, ierr)

  maxRadiusBurned = sqrt(global_maxrsqburned)

  ! if (mesh_mype==MASTER_PE) print*, 'max radius burned material:', maxRadiusBurned


  ! default is to derefine unless a criteria requires resolution
  ! note only mark active blocks, otherwise nothing works
  call Grid_getListOfBlocks(ACTIVE_BLKS, blkList,blkCount)
  do iblk = 1, blkCount
     derefine(blkList(iblk)) = .true.
  end do
  do l = 1,gr_numRefineVars
     iref = gr_refine_var(l)
     ref_cut = gr_refine_cutoff(l)
     deref_cut = gr_derefine_cutoff(l)
     ref_filter = gr_refine_filter(l)
     call gr_markRefineDerefine(mesh_mype,iref,ref_cut,deref_cut,ref_filter)
  end do

#ifdef FLASH_GRID_PARAMESH2
  ! Make sure lrefine_min and lrefine_max are obeyed - KW
  if (gr_numRefineVars .LE. 0) then
     call gr_markRefineDerefine(mesh_mype,-1, 0.0, 0.0, 0.0)
  end if
#endif

  !----------------------------------------------------
  !  2. Impose additional limits on refinement level
  !----------------------------------------------------
  do iblk = 1, blkCount
     ! get info about block
     call Grid_getBlkindexLimits(blkList(iblk), blkLimits, blkLimitsGC)
     call Grid_getBlkPtr(blkList(iblk),solnData,CENTER)

     ! find max enuc and min flff in block + guardcells
     maxenuc = 0.0
     minflff = 1.0
     do k = blkLimitsGC(LOW,KAXIS), blkLimitsGC(HIGH,KAXIS)
        do j = blkLimitsGC(LOW,JAXIS), blkLimitsGC(HIGH,JAXIS)
           do i = blkLimitsGC(LOW,IAXIS), blkLimitsGC(HIGH,IAXIS)
              maxenuc = max(abs(solnData(ENUC_VAR,i,j,k)),maxenuc)
              minflff  = min(solnData(FLFF_MSCALAR,i,j,k),minflff)
           enddo
        enddo
     enddo

     call Grid_releaseBlkPtr(blkList(iblk),solnData)

     !----------------------------------------------------
     !  2.1. impose specified max refinement in fluff
     !     where fluff is determined by flff > FluffThresh (typically 0.99)
     !  regions with minflff > threshold (typcially 0.99) are derefined
     !----------------------------------------------------

     ! check min flff in block (including guardcells)
     !  i.e. block must contain only fluff in order to be prevented from refining
     ! above threshold force derefinement to sim_refFluffLevel
     if ( minflff > sim_refFluffThresh ) then
        if ( lrefine(blkList(iblk)) > sim_refFluffLevel ) then
           refine(blkList(iblk)) = .false.
           derefine(blkList(iblk)) = .true.
        else if ( lrefine(blkList(iblk)) == sim_refFluffLevel ) then
           refine(blkList(iblk)) = .false.
        endif
     ! below threshold but within margin prevent refinement past FluffLevel
     ! but don't force derefinement
     else if ( minflff > sim_refFluffThresh*(1.0-sim_refFluffMargin) &
               .and. lrefine(blkList(iblk)) >= sim_refFluffLevel ) then
        refine(blkList(iblk)) = .false.
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
        if ( lrefine(blkList(iblk)) > sim_refNogenLevel ) then
           refine(blkList(iblk)) = .false.
           derefine(blkList(iblk)) = .true.
        else if ( lrefine(blkList(iblk)) == sim_refNogenLevel ) then
           refine(blkList(iblk)) = .false.
        endif
     ! below threshold but within margin prevent refinement past NogenLevel
     ! but don't force derefinement
     else if ( maxenuc < sim_refNogenEnucThresh &
               .and. lrefine(blkList(iblk)) >= sim_refNogenLevel ) then
        refine(blkList(iblk)) = .false.
     endif
     ! else allow refinement (lrefine_max is enforced below)

     !----------------------------------------------------
     !  2.3. standard refinement min and max
     !----------------------------------------------------

     if (lrefine(blkList(iblk)) < lrefine_min) then
        refine(blkList(iblk)) = .true.
        derefine(blkList(iblk)) = .false.
     else if (lrefine(blkList(iblk)) == lrefine_min) then
        derefine(blkList(iblk)) = .false.
     endif
     if (lrefine(blkList(iblk)) > lrefine_max) then
        refine(blkList(iblk)) = .false.
        derefine(blkList(iblk)) = .true.
     else if (lrefine(blkList(iblk)) == lrefine_max) then
        refine(blkList(iblk)) = .false.
     endif
     
     !----------------------------------------------------
     !  2.4. ejecta phase resolution limits
     !----------------------------------------------------

     if (simTime > sim_refEjectaPhaseStartTime) then

        ejecta_min_cellsize = maxRadiusBurned / sim_refEjectaPhaseMaxRes
        call Grid_getDeltas( blkList(iblk), deltas )
        cellsize = deltas(1)

        ! if cell size for this block is less than min, force derefinement
        if ( cellsize < ejecta_min_cellsize ) then
           derefine(blkList(iblk)) = .true.
           refine(blkList(iblk)) = .false.
        ! only allow refinement if new cell size would be 110% or larger than minimum
        else if ( 0.5*cellsize < 1.1*ejecta_min_cellsize ) then
              refine(blkList(iblk)) = .false.
        endif

     endif

  enddo

  !---------------------------------------------------------------
  ! 3. Do parent-child consistency
  !
  ! for children that are marked derefine, check if parent is 
  ! marked refine and if so unmark derefine
  !---------------------------------------------------------------
  call Grid_getListOfBlocks(ALL_BLKS, blkList,blkCount)
  refine_parent(:) = .false.
  ! open (async) message recieve if parent is off-procssor
  ! otherwise fill directly
  !    message id is child block number on local processor
  nrecv = 0
  do iblk = 1, blkCount
     i = blkList(iblk)
     if (parent(1,i) > 0) then
        if (parent(2,i)/=mesh_mype) then
           nrecv = nrecv+1
           call MPI_IRecv(refine_parent(i), 1, MPI_LOGICAL, parent(2,i), &
                                    i, meshcomm, recvreq(nrecv), ierr)
        else
           refine_parent(i) = refine(parent(1,i))
        endif
     endif
  end do
  ! parents send refine flag to each off-processor child
  nsend = 0
  do iblk = 1, blkCount
     i = blkList(iblk)
     do j = 1,nchild
        if (child(1,j,i) > 0) then
           if (child(2,j,i) /= mesh_mype) then
              nsend = nsend + 1
              call MPI_ISend(refine(i), 1, MPI_LOGICAL, child(2,j,i), &
                                  child(1,j,i), meshcomm, sendreq(nsend), ierr)
           endif
        endif
     enddo
  enddo
  ! wait to recieve all parent info
  if (nrecv > 0) then
     call MPI_Waitall(nrecv,recvreq,recvstat,ierr)
  endif
  ! now reconcile, deferring to parent marked for refine
  do iblk = 1, blkCount
     i = blkList(iblk)
     if (nodetype(i) == LEAF .and. derefine(i) .and. refine_parent(i) ) then
        derefine(i) = .false.
     endif
  enddo
  ! wait until last to ask to have refine() buffer back under our control since
  ! we are not modifying it
  if (nsend > 0) then
     call MPI_Waitall(nsend,sendreq,sendstat,ierr)
  endif


  return
end subroutine Grid_markRefineDerefine
