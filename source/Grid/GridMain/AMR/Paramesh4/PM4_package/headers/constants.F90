!
!!  Licensed under the Apache License, Version 2.0 (the "License");
!!  you may not use this file except in compliance with the License.
!! 
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
! Modification history:
!     Michael L. Rilee, November 2002, *dbz*
!        Initial support for divergenceless prolongation
!     Michael L. Rilee, December 2002, *clean_divb*
!        Support for projecting field onto divergenceless field
!

      module constants

      implicit none
  
      real(kind=kind(1.0d0)), parameter ::  & 
     &        half=0.5d0,  & 
     &        zero=0.0d0, one=1.0d0, two=2.0d0, three=3.0d0, & 
     &        pi=3.141592653589793d0

      end module constants
