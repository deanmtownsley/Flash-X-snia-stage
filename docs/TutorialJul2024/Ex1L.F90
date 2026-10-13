!! Display 1
!! This is how we iterate over blocks in Flash-X.
!! "Tiles" can represent either full blocks (the normal case)
!! or partial blocks ("proper tiles").

  use Grid_iterator, ONLY : Grid_iterator_t

  type(Grid_iterator_t) :: itor

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     ! Do something with the current tile !
     
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

!! The iterator object "itor" is a variable of derived type that contains
!! the current state of the iteration.
!! The details vary highly between different
!! Grid implementations (Paramesh vs Amrex).
!!
!! WHY LIKE THIS? Some reasons:
!! 1. Common interface for Paramesh, Amrex, UG grids
!! 2. Proper tiles can be supported without (much) change to user code
!! 3. Communications (guard cell filling) can be hidden in the intor%next() call

!! Display 2
!! This is how we access the current tile during an interation:


  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)

     ! Do something with the tile with this descriptor !
     
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

!! The tile descriptor "tileDesc" is a variable of derived type that contains
!! information identifying and describing a tile. User code should not be
!! concerend with its internal details. Those details vary between different
!! Grid implementations (Paramesh vs Amrex vs UG).

!! Display 3
!! With the tile descriptor we can access the solution data of a block,
!! using a Fortran POINTER:

  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real, pointer:: solnData(:,:,:,:)

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Do something with the data in solnData !

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

!!  Note that all access to solution data should use pointers, as shown.
!!  Direct access to the underlying data (e.g., Paramesh UNK array) is not supported.
!!  Flash-X does not provide interfaces (like old Grid_getBlkData, Grid_putBlkData)
!!  for copying data to/from user arrays).

!! Display 4
!! For example, this is how solution data are initialized:


  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real, pointer:: solnData(:,:,:,:)

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)

     ! Call the user function to fill in initial data:
     call Simulation_initBlock(solnData, tileDesc)

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

!! Note the different calling convention for Simulation_initBlock!
!! (Passing both solnData and tileDesc is redundant, but is done for convenience.)

!! Display 5
!! We can get metainformation and operate on block data directly:


  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,     ONLY : Grid_tile_t

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)     :: tileDesc
  real, pointer:: solnData(:,:,:,:)

  call Grid_getTileIterator(itor, nodetype=LEAF)
  do while(itor%isValid())
     call itor%currentTile(tileDesc)
     call tileDesc%getDataPtr(solnData, CENTER)

     do k=tileDesc%limits(LOW,KAXIS),tileDesc%limits(HIGH,KAXIS)
        do j=tileDesc%limits(LOW,JAXIS),tileDesc%limits(HIGH,JAXIS)
           do i=tileDesc%limits(LOW,IAXIS),tileDesc%limits(HIGH,IAXIS)
              solnData(DENS_VAR,i,j,k) = 0.001
           end do
        end do
     end do

     call tileDesc%releaseDataPtr(solnData, CENTER)
     call itor%next()
  end do
  call Grid_releaseTileIterator(itor)

!! The solnData array always includes room for guard cells and interior cells.
!!
!! WARNING, GLOBAL INDEXING!
!! The valid index ranges of solnData are different from FLASH!
!! A global range for the indexes (at a given refinement level) is used.
!! The data for each block are if declared by, for example:
!!    real,dimension(NUNK_VARS,1-4:8+4,1-4:8+4,1) :: solnData ! for the first block,
!!    real,dimension(NUNK_VARS,9-4:16+4,1-4:8+4,1) :: solnData ! for the first block,
!!    ...
!! for a 2D case with NXB=NYB=8, NGUARD=4 .
!!
!! WARNING, INDEX REORDERING!
!! This is optional with Paramesh (see setup command), but REQUIRED with Amrex.
!! The actual index order may be converted from
!!    solnData(nvar,i,j,k)
!! in user code to
!!    solnData(i,j,k,nvar)
!! when the source file contains a directive like
!!  !!REORDER(4):solnData
!!





!!!! LEFTOVER !!!!
do k=limits(LOW,KAXIS),limits(HIGH,KAXIS)
   do j=limits(LOW,JAXIS),limits(HIGH,JAXIS)
      do i=limits(LOW,IAXIS),limits(HIGH,IAXIS)
         solnData(DENS_VAR,i,j,k) = 0.001
      end do
   end do
end do






gg


hh




g

