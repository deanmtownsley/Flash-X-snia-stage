!!****ih* source/physics/sourceTerms/Burn/BurnIntegrate/bnIntegrate_interface
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
!!  SYNOPSIS
!!
!!   use bnIntegrate_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the BurnIntegrate module that defines its
!! public interfaces.
!!
!!***
Module bnIntegrate_interface 

   use bnNetwork_interface, ONLY: derivs_t, jakob_t, bjakob_t, steper_t, steper_state_t
  !! These routines are found in the physics/sourceTerms/Burn/BurnIntegrate
  !! They are all used for integrating the ODE of nuclear burning

  interface
     subroutine bn_netIntegrate(btemp,start,stptry,stpmin,stopp,bc, &
                                eps,dxsav,kmax,                     &
                                xrk,yrk,xphys,yphys,xlogi,ylogi,    &
                                nok,nbad,kount,odescal,iprint,      &
                                nrat, ratdum,                       &
                                derivs,jakob,bjakob,steper)
       import :: derivs_t, jakob_t, bjakob_t, steper_t
       implicit none
       integer, intent(IN)  :: xphys,yphys,xlogi,ylogi,nrat
       integer, intent(IN)  :: kmax, iprint
       real, intent(IN)     :: odescal, dxsav, eps
       real, intent(IN)     :: btemp, start, stptry, stpmin, stopp
       real, intent(IN), dimension(nrat)     :: ratdum
       real, intent(INOUT), dimension(yphys) :: bc
       integer, intent(OUT) :: nok, nbad, kount
       real, intent(OUT), dimension(xphys)       :: xrk
       real, intent(OUT), dimension(yphys,xphys) :: yrk
       procedure(derivs_t) :: derivs
       procedure(jakob_t) :: jakob
       procedure(bjakob_t) :: bjakob
       procedure(steper_t) :: steper

       !! bn_netIntegrate calls four external functions derivs, jakob, bjakob, steper
       !! Their interfaces are given below.  They should correspond to (in the case of Aprox13)
       !! derivs = bn_network = aprox13  -- note these were the approximate names in flash2
       !! jakob  = bn_        = saprox13
       !! bjakob = bn_        = baprox13
       !! steper = bn_baderMa28, for example

     end subroutine bn_netIntegrate

  end interface

  interface
     subroutine bn_baderStepMa28(btemp,y,dydx,dfdy,ratdum,nmax,n,xs,htot,nnstep,yout,  & 
          &                  nzo,a,naij,ivect,jvect,jloc,ikeep,iw,w,flag,  & 
          &                  derivs) 
       import :: derivs_t
       implicit none
       integer, intent(IN) :: n, naij, nmax
       integer, intent(INOUT) :: ikeep(1), jloc(naij)  ! because LBR can't figure out what it does
       integer, intent(INOUT) :: iw(1), flag !used by ma28 LU decomposition
       real, intent(OUT)   :: a(naij)
       real, intent(IN)    :: y(n), dydx(n), dfdy(naij), ratdum(:), xs, htot, btemp
       integer, intent(IN) :: nnstep, nzo, ivect(naij), jvect(naij)
       real, intent(INOUT) :: yout(n), w(1)
       procedure(derivs_t) :: derivs
     end subroutine bn_baderStepMa28
  end interface

  interface
     subroutine bn_baderStepGift(btemp,y,dydx,dfdy,ratdum,nmax,n,xs,htot,nnstep,yout,  & 
          &                      derivs)
       import :: derivs_t
       implicit none
       integer, intent(IN) :: n, nmax
       integer, intent(IN) :: nnstep
       real, intent(INOUT) :: yout(n)
       real, intent(IN)    :: btemp, xs, htot, y(n), dydx(n), dfdy(nmax,nmax), ratdum(:)
       procedure(derivs_t) :: derivs
     end subroutine bn_baderStepGift
  end interface

  ! state for bn_pzExtr
  type :: pzExtr_state_t
    real :: x(13)
    real :: qcol(50, 13)
    real :: d(50)
  end type pzExtr_state_t

  interface
     subroutine bn_pzExtr(state,iest,xest,yest,yz,dy,nv) 
       import :: pzExtr_state_t
       implicit none
       type(pzExtr_state_t), intent(IN OUT) :: state
       integer, intent(IN)  :: nv, iest
       real, intent(IN)     :: xest
       real, intent(IN), dimension(nv) :: yest
       real, intent(OUT), dimension(nv) :: dy, yz
     end subroutine bn_pzExtr
  end interface

!-------------------------------------------------------------------------
  interface
     subroutine bn_baderMa28(state,y,dydx,ratdum,nv,x,btemp,htry,eps,yscal,hdid,hnext, & 
          &                   derivs,jakob,bjakob)
       import :: derivs_t, jakob_t, bjakob_t, steper_state_t
       implicit none
       type(steper_state_t), intent(IN OUT) :: state
       integer, intent(IN) :: nv
       real, intent(IN)    :: dydx(nv), yscal(nv), ratdum(:), htry, eps, btemp
       real, intent(INOUT) :: x, y(nv)
       real, intent(OUT)   :: hdid, hnext
       procedure(derivs_t) :: derivs
       procedure(jakob_t) :: jakob
       procedure(bjakob_t) :: bjakob
     end subroutine bn_baderMa28
  end interface

  interface
     subroutine bn_baderGift(state,y,dydx,ratdum,nv,x,btemp,htry,eps,yscal,hdid,hnext, & 
          &                       derivs,jakob,bjakob)
       import :: derivs_t, jakob_t, bjakob_t, steper_state_t
       implicit none
       type(steper_state_t), intent(IN OUT) :: state
       integer, intent(IN) :: nv
       real, intent(IN)    :: dydx(nv), yscal(nv), ratdum(:), htry, eps, btemp
       real, intent(INOUT) :: x, y(nv)
       real, intent(OUT)   :: hdid, hnext
       procedure(derivs_t) :: derivs
       procedure(jakob_t) :: jakob
       procedure(bjakob_t) :: bjakob
     end subroutine bn_baderGift
  end interface


  interface
     subroutine bn_rosenMa28(state,y,dydx,ratdum,n,x,btemp,htry,eps,yscal,hdid,hnext,  & 
          &                      derivs,jakob,bjakob) 
       import :: derivs_t, jakob_t, bjakob_t, steper_state_t
       implicit none
       type(steper_state_t), intent(IN OUT) :: state
       integer, intent(IN) :: n
       real, intent(IN)    :: dydx(n), yscal(n), ratdum(:), htry, eps, btemp
       real, intent(INOUT) :: x, y(n)
       real, intent(OUT)   :: hdid, hnext
       procedure(derivs_t) :: derivs
       procedure(jakob_t) :: jakob
       procedure(bjakob_t) :: bjakob
     end subroutine bn_rosenMa28
  end interface

  interface
     subroutine bn_rosenGift(state,y,dydx,ratdum,n,x,btemp,htry,eps,yscal,hdid,hnext,  & 
          &                      derivs,jakob,bjakob) 
       import :: derivs_t, jakob_t, bjakob_t, steper_state_t
       implicit none
       type(steper_state_t), intent(IN OUT) :: state
       integer, intent(IN) :: n
       real, intent(IN)    :: dydx(n), yscal(n), ratdum(:), htry, eps, btemp
       real, intent(INOUT) :: x, y(n)
       real, intent(OUT)   :: hdid, hnext
       procedure(derivs_t) :: derivs
       procedure(jakob_t) :: jakob
       procedure(bjakob_t) :: bjakob
     end subroutine bn_rosenGift
  end interface



end Module bnIntegrate_interface
