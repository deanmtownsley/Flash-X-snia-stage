!!****ih* source/Grid/GridMain/AMR/Paramesh4/block_metadata
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!!
!!
!!
!!****

#include "constants.h"
#include "Simulation.h"

module block_metadata

    implicit none

    private

    type, public :: block_metadata_t
        integer :: id
        integer :: cid(MDIM)
        integer :: stride(MDIM)
        integer :: level
        integer :: limits(LOW:HIGH, MDIM)
        integer :: limitsGC(LOW:HIGH, MDIM)
        integer :: localLimits(LOW:HIGH, MDIM)
        integer :: localLimitsGC(LOW:HIGH, MDIM)
      contains
        procedure,public :: getCellCoords
    end type block_metadata_t

    contains
      subroutine getCellCoords(self, coords, axis, edge, indexing, offset, npts)
        use Grid_data, ONLY : gr_globalDomain, gr_delta, gr_lrefineMax


        class(block_metadata_t),intent(IN) :: self
        real,intent(OUT) :: coords(:)
        integer,intent(IN) :: axis
        integer,intent(IN), OPTIONAL :: edge
        integer,intent(IN), OPTIONAL :: indexing
        integer,intent(IN), OPTIONAL :: npts
        integer,intent(IN), OPTIONAL :: offset

        integer :: edgeEff
        integer :: indexingEff
        integer :: nptsEff
        integer :: offsetEff
        integer :: i, first

        if (present(edge)) then
           edgeEff = edge
        else
           edgeEff = CENTER
        end if
        if (present(indexing)) then
           indexingEff = indexing
        else
           indexingEff = GLOBALIDX1
        end if
        if (present(offset)) then
           offsetEff = offset
        else
           offsetEff = 0
        end if
        if (present(npts)) then
           nptsEff = npts
        else
           nptsEff = size(coords)
        end if

        select case (indexingEff)
        case(EXTERIOR)
           first=self%cid(axis)-1-self%stride(axis)*NGUARD
        case(INTERIOR)
           first=self%cid(axis)-1
        case default
           first=0
        end select

        first = first + offsetEff*self%stride(axis)
        if((edge==CENTER).and.(self%stride(axis)==1))then
           do i = 1,nptsEff
              coords(i)= gr_globalDomain(LOW,axis) + (first+0.5)*gr_delta(axis,gr_lrefineMax)
              first=first+self%stride(axis)
           end do
        else
           if(edge==RIGHT_EDGE)first=first+self%stride(axis)
           if((edge==CENTER).and.(self%stride(axis)>1))first=first+self%stride(axis)/2
           do i = 1,nptsEff
              coords(i)= gr_globalDomain(LOW,axis) + first*gr_delta(axis,gr_lrefineMax)
              first=first+self%stride(axis)
           end do
        end if

      end subroutine getCellCoords
end module block_metadata

