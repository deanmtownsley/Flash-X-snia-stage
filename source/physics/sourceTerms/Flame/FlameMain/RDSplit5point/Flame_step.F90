!!****if* source/physics/sourceTerms/Flame/FlameMain/RDSplit5point/Flame_step
!!
!! NAME
!!
!!  Flame_step
!!
!! SYNOPSIS
!!
!!  call Flame_step( real(in)    :: dt)
!!
!! DESCRIPTION
!!
!!  see Flame_interface.F90 for function description
!!
!! ARGUMENTS
!!
!!   num_blocks : 
!!
!!   blocklist : 
!!
!!   dt : 
!!
!! AUTOGENROBODOC
!!
!!
!!***

! see Flame_interface.F90 at top level for function description
!
! Dean Townsley 2008
!

! Implementation details
!
! This implementation is operator split betwee reaction and diffusion
! but dimensionally unsplit.  Diffusion is treated by directly computing
! the Laplacian instead of differencing fluxes, so that this is not a
! conservative diffusion operator.
!
! Computation of the Laplacian itself is done in a subroutine, as its
! form depends heavily on the mesh geometry

!!REORDER(4): solnData

#define DEBUG_GRID_GCMASK

#include "Simulation.h"
#include "constants.h"
subroutine Flame_step( dt )    

  use Flame_data

  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,      ONLY : Grid_tile_t
  use Grid_interface, ONLY : Grid_getTileIterator,Grid_fillGuardCells
  use fl_fsInterface, only : fl_flameSpeed
  use fl_effInterface, only: fl_effects
  use fl_interface, only : fl_laplacian
  use Driver_interface, only : Driver_abort
  use Timers_interface, only : Timers_start, Timers_stop
  use Logfile_interface, only : Logfile_stamp, Logfile_stampVarMask
       
  implicit none
  real,    INTENT(in)                        :: dt

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  integer, dimension(2,MDIM) :: fspeedLimits
  integer :: istat
  integer :: n, bid

  real, pointer, dimension(:,:,:,:) :: solnData
  real, allocatable, dimension(:,:,:) :: flam, flamdot, flamespeed, lapl

  real :: f, inv_dt
  integer :: i,j,k
  integer :: sizeI, sizeJ, sizeK

  if( .not. fl_useFlame ) return

  call Timers_start("flame")

  inv_dt = 1.0/dt

  if (fl_gcDoLogMask) then
     call Logfile_stamp('calling guardcell fill with mask logging on','[Flame_step]')
#ifdef DEBUG_GRID_GCMASK
     call Logfile_stampVarMask(fl_gcMask, fl_gcDoEos, '[Flame_step]', 'gcMask')
#endif
  end if
  call Grid_fillGuardCells(CENTER, ALLDIR, eosMode=MODE_DENS_EI, &
    doEos=fl_gcDoEos, maskSize=fl_gcMaskSize, mask=fl_gcMask, &
    makeMaskConsistent=.false., doLogMask=fl_gcDoLogMask)
  fl_gcDoLogMask=.false.

  call Grid_getTileIterator( itor, nodetype=LEAF )

  do while( itor%isValid() )
     call itor%currentTile(tileDesc)

     call tileDesc%getDataPtr( solnData, CENTER)

     allocate( flam( tileDesc%limits(LOW,IAXIS) : tileDesc%limits(HIGH,IAXIS), &
                     tileDesc%limits(LOW,JAXIS) : tileDesc%limits(HIGH,JAXIS), &
                     tileDesc%limits(LOW,KAXIS) : tileDesc%limits(HIGH,KAXIS)), &
               STAT=istat )
     if (istat /= 0) call Driver_abort("Cannot allocate flam in Flame_step")
     allocate( flamdot( tileDesc%limits(LOW,IAXIS) : tileDesc%limits(HIGH,IAXIS), &
                        tileDesc%limits(LOW,JAXIS) : tileDesc%limits(HIGH,JAXIS), &
                        tileDesc%limits(LOW,KAXIS) : tileDesc%limits(HIGH,KAXIS)), &
               STAT=istat )
     if (istat /= 0) call Driver_abort("Cannot allocate flamdot in Flame_step")
     allocate( flamespeed( tileDesc%limits(LOW,IAXIS) : tileDesc%limits(HIGH,IAXIS), &
                           tileDesc%limits(LOW,JAXIS) : tileDesc%limits(HIGH,JAXIS), &
                           tileDesc%limits(LOW,KAXIS) : tileDesc%limits(HIGH,KAXIS)), &
               STAT=istat )
     if (istat /= 0) call Driver_abort("Cannot allocate flamespeed in Flame_step")
     allocate( lapl( tileDesc%limits(LOW,IAXIS) : tileDesc%limits(HIGH,IAXIS), &
                     tileDesc%limits(LOW,JAXIS) : tileDesc%limits(HIGH,JAXIS), &
                     tileDesc%limits(LOW,KAXIS) : tileDesc%limits(HIGH,KAXIS)), &
               STAT=istat )
     if (istat /= 0) call Driver_abort("Cannot allocate lapl in Flame_step")

     ! extract flam variable
     ! need two layers in GCs becausee of RD splitting
     flam( tileDesc%limits(LOW,IAXIS)-2 : tileDesc%limits(HIGH,IAXIS)+2 , &
           tileDesc%limits(LOW,JAXIS)-2*K2D : tileDesc%limits(HIGH,JAXIS)+2*K2D , &
           tileDesc%limits(LOW,KAXIS)-2*K3D : tileDesc%limits(HIGH,KAXIS)+2*K3D) &
         = solnData(FLAM_MSCALAR, tileDesc%limits(LOW,IAXIS)-2 : tileDesc%limits(HIGH,IAXIS)+2 , &
                 tileDesc%limits(LOW,JAXIS)-2*K2D : tileDesc%limits(HIGH,JAXIS)+2*K2D , &
                 tileDesc%limits(LOW,KAXIS)-2*K3D : tileDesc%limits(HIGH,KAXIS)+2*K3D)

     call fl_flameSpeed(solnData, flamespeed, tileDesc, 2)

     do k = tileDesc%limits(LOW,KAXIS)-2*K3D, tileDesc%limits(HIGH,KAXIS)+2*K3D
        do j = tileDesc%limits(LOW,JAXIS)-2*K2D, tileDesc%limits(HIGH,JAXIS)+2*K2D
           do i = tileDesc%limits(LOW,IAXIS)-2, tileDesc%limits(HIGH,IAXIS)+2
              f = flam(i,j,k)
              flam(i,j,k) = f + dt*fl_R_over_s*flamespeed(i,j,k)*(f-fl_epsilon_0)*(1.0+fl_epsilon_1-f)
           enddo
        enddo
     enddo

     ! 1 specifies the step size should be 1 grid cell
     ! cannot be any larger because flam is filled with only 2 guard cell layers
     call fl_laplacian(lapl, flam, 1, tileDesc)

     do k = tileDesc%limits(LOW,KAXIS), tileDesc%limits(HIGH,KAXIS)
        do j = tileDesc%limits(LOW,JAXIS), tileDesc%limits(HIGH,JAXIS)
           do i = tileDesc%limits(LOW,IAXIS), tileDesc%limits(HIGH,IAXIS)
              flam(i,j,k) = max(0.0, min(1.0, flam(i,j,k) + dt*fl_kappa_over_s*flamespeed(i,j,k)*lapl(i,j,k) ) )
              flamdot(i,j,k) = (flam(i,j,k) - solnData(FLAM_MSCALAR,i,j,k))*inv_dt
              solnData(FLAM_MSCALAR, i,j,k) = flam(i,j,k)
           enddo
        enddo
     enddo

     deallocate(lapl)
     deallocate(flamespeed)
     deallocate(flam)

     call fl_effects( solnData, flamdot, dt, tileDesc)

     deallocate(flamdot)

     call tileDesc%releaseDataPtr(solnData, CENTER)

  enddo

  call Timers_stop("flame")

  return
end subroutine Flame_step
