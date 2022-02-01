!!****if* source/IO/IOParticles/hdf5/serial/io_ptWriteParticleData
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
!! io_ptWriteParticleData
!!
!!
!! SYNOPSIS
!!
!! call io_ptWriteParticleData(integer(io_fileID_t)(in) :: fileID,
!!                     integer(in) :: globalNumParticles,
!!                     integer(in) :: localNumParticles,
!!                     integer(in) :: particleOffset,
!!                     character(in)(len=OUTPUT_PROP_LENGTH):: partAttributeLabels(NPART_PROPS)
!!                     logical(in) :: particlesToCheckpoint))
!!
!!
!!
!! DESCRIPTION
!!
!!    This routine writes out the particle data in a separate hdf5 file
!!    It calls  ioh5_write_particles
!!    This is being done to make particles easy to debug and since the particles
!!    are not associated with the mesh data it makes since to separate it from
!!    the checkpoint files and plotfiles
!!
!!
!! ARGUMENTS 
!!
!!   
!!   fileID - file handle to an open file, in this case an HDF5 file identifier
!!   
!!   globalNumParticles - the total number of particles in the simulation
!!   
!!   localNumParticles - number of particles on a given processor
!!   
!!   particleOffset - global particle offset, in this IO implementation
!!                    particleOffset is not used.  Each proc writes its own
!!                    particles to a unique file (either checkpoint or particle plotfile)
!!                    so no global offset is needed.  (Kept in the interface for
!!                    consistency with other IO implementations like pnetcdf and hdf5.)
!!
!!   partAttributeLabels - string label names for the particle properties.   
!!                         (like 'xpos', 'yvel', 'tag' etc)
!!
!!   particlesToCheckpoint - logical value indicating if particles should be written
!!                           to a checkpoint file or to a plotfile. (not used in this 
!!                           implementation)
!!
!!
!! NOTES
!!   To control particle output there are a few different parameters
!!   
!!
!!
!!***

#include "Simulation.h"

subroutine io_ptWriteParticleData( fileID, globalNumParticles, &
     localNumParticles, particleOffset, partAttributeLabels, particlesToCheckpoint)
  
  
  use IO_data, ONLY : io_realParmNames, io_realParmValues, io_numRealParms, &
       io_intParmNames, io_intParmValues, io_numIntParms, &
       io_logParmNames, io_logParmValues, io_numLogParms, &
       io_strParmNames, io_strParmValues, io_numStrParms, &
       io_realScalarNames, io_realScalarValues, io_numRealScalars, &
       io_intScalarNames, io_intScalarValues, io_numIntScalars, &
       io_logScalarNames, io_logScalarValues, io_numLogScalars, &
       io_strScalarNames, io_strScalarValues, io_numStrScalars, &
       io_logToIntScalarValues, io_logToIntParmValues, &
       io_geometry, io_setupCall, io_buildDir, io_flashRelease, &
       io_fileCreationTime, io_buildDate, io_buildMachine, io_cflags, &
       io_fflags, io_setupTimeStamp, io_buildTimeStamp, io_outputSplitNum, &
       io_fileFormatVersion, io_globalMe, io_globalNumProcs

  use io_intfTypesModule, ONLY : io_fileID_t
  use Grid_interface, ONLY : Grid_getLocalNumBlks, Grid_sortParticles, Grid_getTileIterator, Grid_releaseTileIterator
  use Particles_data, ONLY : particles, pt_maxPerProc
  use Grid_data, ONLY : gr_globalNumBlocks

#ifdef FLASH_GRID_AMREX
  use Grid_interface, ONLY : Grid_countParticles, Grid_countParticlesByBlock
  use Particles_interface, ONLY : Particles_copyFromMeshOwned
  use Particles_data, ONLY : pt_containers
#endif
  
  use Grid_iterator, ONLY : Grid_iterator_t
  use Grid_tile,        ONLY : Grid_tile_t

  implicit none     

#include "constants.h"
#include "Flashx_mpi.h"

  integer(io_fileID_t), intent(in) ::  fileID
  integer, intent(in) :: globalNumParticles, localNumParticles, particleOffset
  character (len=OUTPUT_PROP_LENGTH), intent(in) :: partAttributeLabels(NPART_PROPS)
  logical, intent(in) :: particlesToCheckpoint

 ! for logfile output
  character(len=MAX_STRING_LENGTH), dimension(2,2) :: strBuff
  integer :: ierr, blocksPerFile, localNumBlocks, pstart, pend, count
  real, allocatable :: particlest(:,:)

  integer :: jproc, localnumparticlest, status(MPI_STATUS_SIZE)
  integer :: blockOffset, localNumBlockst, lb, beginIndex, endIndex
  integer,parameter :: particleTypes=NPART_TYPES
  integer,dimension(MAXBLOCKS,particleTypes) :: particlesPerBlkt, particlesPerBlk
  integer :: partOffset
  integer :: l_numParticles, particleType

  integer :: p_count, ind, t_blk

  type(Grid_iterator_t) :: itor
  type(Grid_tile_t)    :: tileDesc


  l_numParticles=localNumParticles


  allocate(particlest(NPART_PROPS, pt_maxPerProc))

  call Grid_getLocalNumBlks(localNumBlocks)

!! This call returns particles sorted by block numbers and also the 
!! the number of particles that each block has.

   ! set to zero
   particlesPerBlk(:,:) = 0.0
   particlesPerBlkt(:,:) = 0.0

#ifdef FLASH_GRID_AMREX
  call Grid_countParticles(NPART_PROPS,l_numParticles,particleTypes, pt_maxPerProc)
! set ParticlesPerBlk for AMReX
  call Grid_countParticlesByBlock(particlesPerBlk)
#else
  call Grid_sortParticles(particles,NPART_PROPS,l_numParticles, particleTypes,&
       pt_maxPerProc,particlesPerBlk, BLK_PART_PROP)  
#endif

  if((io_globalMe == MASTER_PE) .and. (.not. particlesToCheckpoint)) then
    
    call io_h5write_header(io_globalMe, NPART_PROPS, fileID, io_geometry, &
         partAttributeLabels, io_setupCall, io_fileCreationTime, &
         io_flashRelease, io_buildDate, io_buildDir, io_buildMachine, &
         io_cflags, io_fflags, io_setupTimeStamp, io_buildTimeStamp, &
         io_fileFormatVersion, io_outputSplitNum)

    call io_prepareListsWrite()

         !! write the runtime parameters
     call io_h5write_runtime_parameters(io_globalMe, &
          fileID, &
          io_numRealParms, &
          io_realParmNames, &
          io_realParmValues, &
          io_numIntParms, &
          io_intParmNames, &
          io_intParmValues, &
          io_numStrParms, &
          io_strParmNames, &
          io_strParmValues, &
          io_numLogParms, &
          io_logParmNames, &
          io_logToIntParmValues, &
          io_outputSplitNum)
     
     !! write the scalars
     call io_h5write_scalars(io_globalMe, &
          fileID, &
          io_numRealScalars, &
          io_realScalarNames, &
          io_realScalarValues, &
          io_numIntScalars, &
          io_intScalarNames, &
          io_intScalarValues, &
          io_numStrScalars, &
          io_strScalarNames, &
          io_strScalarValues, &
          io_numLogScalars, &
          io_logScalarNames, &
          io_logToIntScalarValues, &
          io_outputSplitNum)

     call io_finalizeListsWrite()

  end if !(io_globalMe == MASTER_PE) .and. (.not. particlesToCheckpoint)

  !-----------------------------------------------------------------------------
  ! loop over all of the processors.  All the data is moved to processor 0 for
  ! storage using MPI sends and receives.
  !-----------------------------------------------------------------------------
  partOffset = 0  !in serial version, better to calculate than to communicate
  blockOffset = 0

  do jproc = 0,io_globalNumProcs-1
     
     if (io_globalMe == MASTER_PE) then
        ! fetch localNumblocks from other processors
        if (jproc /= 0) then
          call MPI_RECV (localNumParticlest,1,FLASH_INTEGER,jproc, & 
                1,MPI_COMM_WORLD,status,ierr)
          if (localNumParticlest > 0) then
             call MPI_RECV(particlest(1,1), NPART_PROPS*localNumParticlest, FLASH_REAL, &
                      jproc, 2, MPI_COMM_WORLD, status, ierr)
           end if !(if localNumParticlest > 0)
           
           call MPI_RECV (localNumBlockst,1,FLASH_INTEGER,jproc, & 
                3,MPI_COMM_WORLD,status,ierr)
           
           if (localNumBlockst > 0) then
              call MPI_RECV (particlesPerBlkt(1,1), localNumBlockst,FLASH_INTEGER,jproc, & 
                   4,MPI_COMM_WORLD,status,ierr)
           end if
      else
          localNumParticlest = localNumParticles

#ifdef FLASH_GRID_AMREX
          count = 0
          call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
           do while(itor%isValid())
               call itor%currentTile(tileDesc)
               do ind=1,NPART_TYPES
                 p_count = pt_containers(ind)%num_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
                 pstart = count + 1
                 pend =  count + p_count
                 particleType = ind
                 call Particles_copyFromMeshOwned(p_count, particleType, & 
                    pt_containers(ind)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index) , &
                                             particlest(:, pstart:pend))
                 count = count + p_count
               enddo
               call itor%next()
            enddo             ! leaf itor enddo

            call Grid_releaseTileIterator(itor)
#else
          particlest(:,1:localNumParticlest) = particles(:,1:localNumParticlest)
#endif          

          localNumBlockst = localNumBlocks
          particlesPerBlkt(1:localNumBlockst,1) = particlesPerBlk(1:localNumBlockst,1)
        end if

        ! Write particles
        if (localNumParticlest > 0) then
           call io_h5write_particles(io_globalMe, &
                fileID, &
                globalNumParticles, &
                localNumParticlest, &
                NPART_PROPS, &
                partOffset, &
                particlest, &
                partAttributeLabels)      
        end if !if localNumParticlest > 0

        if(localNumBlockst > 0) then
           call io_h5write_localnp(io_globalMe, &
                fileID, &
                particlesPerBlkt, &
                localNumBlockst, &
                gr_globalNumBlocks, &
                blockOffset)
        endif
        !--------------------------------------------------------------------
        ! end local block loop
        !--------------------------------------------------------------------
        
        ! increment the global block number -- 
        !we just wrote localNumBlockst blocks from
        ! processor jproc to the output file
        partOffset = partOffset + localNumParticlest
        blockOffset = blockOffset + localNumBlockst

     else ! if (Io_GlobalMe == MASTER_PE)
        if (jproc == Io_GlobalMe) then
           call MPI_SEND(localNumParticles, 1, FLASH_INTEGER, 0, & 
                1, MPI_COMM_WORLD, ierr)
           if (localNumParticles > 0) then
#ifdef FLASH_GRID_AMREX
              count = 0
              call Grid_getTileIterator(itor, LEAF, tiling=.FALSE.)
              do while(itor%isValid())
                 call itor%currentTile(tileDesc)
                 do ind=1,NPART_TYPES
                    p_count = pt_containers(ind)%num_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index)
              
                    pstart = count + 1
                    pend =  count + p_count
                    particleType = ind
  
                    call Particles_copyFromMeshOwned(p_count, particleType, & 
                       pt_containers(ind)%get_particles(tileDesc%level-1,tileDesc%grid_index, tileDesc%tile_index) , &
                       particlest(:, pstart:pend))
                    count = count + p_count
  
                    ! if we knew how to get a count of lost particles separately 
                    ! we'd increment localElementsPerType(1+LOSTBLK..) instead.   
                 enddo
                 call itor%next()
              enddo             ! leaf itor enddo
  
              call Grid_releaseTileIterator(itor)

              call MPI_SEND(particlest(1,1), NPART_PROPS*localNumParticles, FLASH_REAL, &
                   0, 2, MPI_COMM_WORLD, ierr)
#else
              call MPI_SEND(particles(1,1), NPART_PROPS*localNumParticles, FLASH_REAL, &
                   0, 2, MPI_COMM_WORLD, ierr)
#endif
           endif !localNumParticles > 0
           
           call MPI_SEND(localNumBlocks, 1, FLASH_INTEGER, 0, & 
                3, MPI_COMM_WORLD, ierr)

           if (localNumBlocks > 0) then
               call MPI_SEND(particlesPerBlk(1,1), localNumBlocks, FLASH_INTEGER, 0, & 
                    4, MPI_COMM_WORLD, ierr)
           end if

        end if !jproc == Io_GlobalMe

     end if                 ! if Io_GlobalMe == MASTER_PE

     call MPI_BARRIER (MPI_COMM_WORLD, ierr)


     !------------------------------------------------------------------------
     ! end processor loop
     !------------------------------------------------------------------------
  end do
  deallocate(particlest)
  return

end subroutine io_ptWriteParticleData
