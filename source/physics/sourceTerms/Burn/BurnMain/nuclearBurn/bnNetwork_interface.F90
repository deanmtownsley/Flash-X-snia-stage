!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/bnNetwork_interface
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
!!  bnNetwork_interface
!!
!!
!! SYNOPSIS
!!
!!  use bnNetwork_interface, ONLY: 
!!
!!
!! DESCRIPTION
!!
!!  Subroutine argument descriptions
!!
!!***

#include "Simulation.h"

Module bnNetwork_interface

  interface 
     subroutine bn_initNetwork
       implicit none
     end subroutine bn_initNetwork
  end interface

  interface
     subroutine bn_finalizeNetwork
       implicit none
     end subroutine bn_finalizeNetwork
  end interface

  interface
     subroutine bn_networkRates(btemp, bden, abar, zbar, z2bar, ytot1, bye, nrat, ratraw)
        implicit none
        integer, intent(IN) :: nrat
        real, intent(IN) :: btemp, bden
        real, intent(IN) :: abar, zbar, z2bar, ytot1, bye
        real, dimension(nrat), intent(OUT) :: ratraw
     end subroutine bn_networkRates
  end interface

  interface
     subroutine bn_networkTable(btemp, bden, abar, zbar, z2bar, ytot1, bye, &
                                nrat, nrattab, &
                                rattab, ttab, dtab, ratraw)
        implicit none
        integer, intent(IN) :: nrat, nrattab
        real, intent(IN) :: btemp, bden, abar, zbar, z2bar, ytot1, bye
        real, intent(OUT), dimension(nrattab) :: ttab
        real, intent(OUT), dimension(nrat, nrattab) :: rattab
        real, intent(OUT), dimension(nrat) :: dtab
        real, intent(OUT), dimension(nrat) :: ratraw
     end subroutine bn_networkTable
  end interface

  interface 
     subroutine bn_gift(ab,n1,n2)
       implicit none
       integer, INTENT(in) :: n1,n2
       real, INTENT(inout) :: ab(n1,n2)
     end subroutine bn_gift
  end interface

  interface
     subroutine bn_networkScreen(screen_init, btemp, bden, bye, ratraw, ymass, scfac, nrat, ratdum)
        implicit none
        logical, intent(INOUT) :: screen_init
        integer, intent(IN) :: nrat
        real, intent(IN) :: btemp, bden, bye, ymass(NSPECIES)
        real, intent(IN OUT) :: scfac(nrat), ratraw(nrat), ratdum(nrat)
     end subroutine bn_networkScreen
  end interface

  interface
     subroutine bn_networkWeak(y)
       implicit none
       real, intent(IN) :: y(NSPECIES)
     end subroutine bn_networkWeak
  end interface



! abstract interfaces for bn_netIntegrate
  abstract interface
     subroutine derivs_t(tt,y,btemp,ratdum,dydt)
       implicit none
       real, intent(IN)  :: tt, y(:), btemp, ratdum(:)     ! NOTE: btemp is reserved for Aprox19
       real, intent(OUT) :: dydt(:)
     end subroutine derivs_t
   end interface

   abstract interface
     subroutine jakob_t(tt,y,btemp,ratdum,dfdy,nzo,nDummy) ! = bn_networkSparseJakob or bn_networkDenseJakob
       implicit none                                       ! See notes in bnNetwork_interface about this...
       integer, intent(IN) :: nzo, nDummy
       real, intent(IN)    :: tt, btemp, ratdum(:)
       real, intent(INOUT) :: y(:)
       real, intent(OUT)   :: dfdy(:,:)
     end subroutine jakob_t
   end interface

   abstract interface
     subroutine bjakob_t(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  ::   iloc(:),jloc(:),np
       integer, intent(OUT) ::   nzo
     end subroutine bjakob_t
   end interface

   ! for steper_t, we need a state variable to "save" the state between steper calls.
   ! NOTE: this will only work with *Gift steppers
   type steper_state_t
      logical :: ifirst, first, reduct
      integer :: nvold
      real    :: eps1, epsold
      integer :: kmax, kopt
      integer :: nseq(8)
      real    :: a(8)
      real    :: alf(7, 7)
   end type steper_state_t

   abstract interface
     subroutine steper_t(state,y,dydx,ratdum,nv,x,btemp,htry,eps,yscal,hdid,hnext, &
                         derivs,jakob,bjakob)
       import :: steper_state_t, derivs_t, jakob_t, bjakob_t
       implicit none
       type(steper_state_t), intent(IN OUT) :: state
       integer, intent(IN) :: nv
       real, intent(INOUT) :: y(nv)
       real, intent(IN)    :: dydx(nv), yscal(nv), ratdum(:), htry, eps, btemp
       real, intent(OUT)   :: hdid, hnext
       real, intent(INOUT) :: x
       procedure(derivs_t) :: derivs
       procedure(jakob_t) :: jakob
       procedure(bjakob_t) :: bjakob
     end subroutine steper_t
  end interface
!---------------------------------------------------


!----------------------------------------------------
  ! NOTE: should have the same signature as derivs_t
  interface
     subroutine bn_network(tt,y,btemp,ratdum,dydt)
       implicit none
       real, intent(IN)  :: tt, y(:), btemp, ratdum(:)     ! NOTE: btemp is reserved for Aprox19
       real, intent(OUT) :: dydt(:)
     end subroutine bn_network
  end interface


!-----------------------
!  These two routines are supposed to have the same interface.  Sadly, they
!     do not -- dfdy has different ranks.  bn_networkSparseJakob takes
!     the last argument as a fixed index for the second rank of dfdy
!     to treat dfdy as the same rank of bn_networkDenseJakob.
!  Normally, the actual argument for the nDummy in the bn_networkSparseJakob
!     should be an integer parameter with value of one.
!  NOTE: below two subroutines' procedures are provided as jakob_t
  interface
     subroutine bn_networkSparseJakob(tt,y,btemp,ratdum,dfdy,nzo,nDummy)
       implicit none
       integer, intent(IN) :: nzo, nDummy ! dummy index for dealing with dfdy
       real, intent(IN)    :: tt, btemp, ratdum(:)
       real, intent(INOUT) :: y(:)
       real, intent(OUT)   :: dfdy(:,:)   ! NOTE: it only updates dfdy(:,nDummy)
     end subroutine bn_networkSparseJakob
  end interface

  interface 
     subroutine bn_networkDenseJakob(tt,y,btemp,ratdum,dfdy,nlog,nphys)
       implicit none
       integer, intent(IN) :: nlog, nphys
       real, intent(IN)    :: tt, btemp, ratdum(:)
       real, intent(INOUT) ::  y(:)
       real, intent(OUT)   ::  dfdy(:,:)
     end subroutine bn_networkDenseJakob
  end interface


!----------------------------------------------
! This is the routine passed as bjakob. Bizarrely enough, there is no
!!  equivalent of bn_networkDensePointers
  ! NOTE: should have the same signature as bjakob_t
  interface
     subroutine bn_networkSparsePointers(iloc,jloc,nzo,np)
       implicit none
       integer, intent(IN)  ::   iloc(:),jloc(:),np
       integer, intent(OUT) ::   nzo
     end subroutine bn_networkSparsePointers
  end interface


end Module bnNetwork_interface

