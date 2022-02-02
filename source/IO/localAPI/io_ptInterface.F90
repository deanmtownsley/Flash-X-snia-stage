!!****ih* source/IO/localAPI/io_ptInterface
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
!!  io_ptInterface
!!
!! SYNOPSIS
!!  use io_ptInterface
!!
!! DESCRIPTION
!!
!! This is the interface module for particle IO.  This handles whatever
!! correctoions need to be placed in checkpoint files to the next particle 
!! output time, as well as initialzation of particles.  These must only be
!! called from within IO.
!! 
!!***

#include "constants.h"
#include "Simulation.h"

module io_ptInterface

  interface 
     subroutine io_ptCorrectNextPartTime(simTime, oldTime)
        real, intent(IN) :: simTime
        real, intent(out) :: oldTime
     end subroutine io_ptCorrectNextPartTime
  end interface

  interface
     subroutine io_ptResetNextFile(savedNext)
       real, intent(IN) :: savedNext
     end subroutine io_ptResetNextFile
  end interface

  interface
     subroutine io_ptInit()
       
       implicit none 
       
     end subroutine io_ptInit
  end interface

  interface
     subroutine io_ptSendOutputData()
       implicit none

     end subroutine io_ptSendOutputData
  end interface

  interface
     subroutine io_ptWriteParticleData( fileID, &
          globalNumParticles, localNumParticles, particleOffset, &
          partAttributeLabels, particlesToCheckpoint)
       use io_intfTypesModule, ONLY : io_fileID_t
       implicit none
       integer(io_fileID_t),intent(IN) :: fileID
       integer, intent(IN) :: globalNumParticles, &
            localNumParticles, particleOffset
       character (len=OUTPUT_PROP_LENGTH), intent(IN) :: partAttributeLabels(NPART_PROPS)
       logical, intent(IN):: particlesToCheckpoint
     end subroutine io_ptWriteParticleData
  end interface

  interface
     subroutine io_ptReadParticleData()
       implicit none
     end subroutine io_ptReadParticleData
  end interface
       
  interface
     subroutine io_ptCreateSubset ( subsetIndex, numProperties, &
          numParticles, partLabels, particlest, &
          subsetSize, subsetName, subsetLabelName, subsetLabels, &
          moreSubsetsRemain)   
       implicit none
       integer, intent(IN) :: subsetIndex, numProperties, &
            numParticles
       character(len=OUTPUT_PROP_LENGTH), dimension(numProperties), &
            intent(IN) :: partLabels
       real, dimension(numProperties,numParticles), intent(INOUT) :: particlest
       integer, dimension(2), intent(OUT) :: subsetSize
       character(len=MAX_STRING_LENGTH), intent(OUT) :: subsetName, subsetLabelName
       character(len=OUTPUT_PROP_LENGTH), dimension(numProperties), &
            intent(OUT) :: subsetLabels
       logical, intent(OUT) :: moreSubsetsRemain
     end subroutine io_ptCreateSubset
  end interface

end module io_ptInterface
