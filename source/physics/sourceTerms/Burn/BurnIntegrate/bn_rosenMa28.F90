!!****if* source/physics/sourceTerms/Burn/BurnIntegrate/bn_rosenMa28
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
!!   bn_rosenMa28
!!
!! SYNOPSIS
!!  subroutine bn_rosenMa28(real(IN)       ::y(:),
!!                          real(IN)       ::dydx(:),
!!                          real(IN)       ::ratdum(:),
!!                          integer(IN)    ::n,
!!                          real(IN)       ::x,
!!                          real(IN)       ::btemp,
!!                          real(IN)       ::htry,
!!                          real(IN)       ::eps,
!!                          real(IN)       ::yscal(:),
!!                          real(OUT)      ::hdid,
!!                          real(OUT)      ::hnext, 
!!                          procedure(IN)  :: derivs,
!!                          procedure(IN)  :: jakob,
!!                          procedure(IN)  :: bjakob) 
!!
!! DESCRIPTION
!!
!!  fourth order rosenbrock step for integrating stiff ode's with monitoring 
!!  of local truncation error to adjust the stepsize. input are the dependent  
!!  variable y(1:n) and its derivative dydx(1:n), at the starting value of the  
!!  independent variable x. also input are the stepsize to be attempted htry,  
!!  the desired fractional accuracy eps, and the vector yscal(1:n) against  
!!  which the error is scaled. derivs is a routine that computes the righthand  
!!  side of the system of ode's. on output the y and x are replaced by their  
!!  new values, hdid is the stepsize that was actually accomplished and hnext  
!!  is the estimate for the next stepsize. 
!!   
!!  nmax is the maximum number of ode's.  
!!  grow and shrnk are the extremes of the time step growth  
!!  errcon = (grow/safety)**(1/pgrow) handles the case when errmax = 0.0 
!!  
!!  routine modified to assume dfdx = 0.0 for reaction networks
!!  
!!  sparse ma28 algebra version
!!
!! ARGUMENTS
!!
!!   y       - dependent variable, array of size y(1:n)
!!   dydx    - derivative of dependent variable, array of size dydx(1:n)
!!   ratdum  - reaction rate
!!   n       - number of dependent variables
!!   x       - independent variable
!!   btemp   - temperature
!!   htry    - attempted stepsize 
!!   eps     - desired fractional accuracy
!!   yscal   - vector of size yscal(1:n) for scaling error
!!   hdid    - stepsize used
!!   hnext   - estimate for next stepsize
!!   derivs  - procedure(IN) name of the routine that contains the odes
!!   jakob   - procedure(IN) name of the routine that contains the jacobian of the odes
!!   bjakob  - procedure(IN) name of the routine that sets the pointers of the sparse jacobian
!!
!!
!!***
!!---------------------------------------------------------------------------------

subroutine bn_rosenMa28(state,y,dydx,ratdum,n,x,btemp,htry,eps,yscal,hdid,hnext,  & 
     &                      derivs,jakob,bjakob) 

  use Driver_interface, ONLY : Driver_abort
  use Logfile_interface, ONLY: Logfile_stampMessage
  use bnIntegrate_interface, ONLY: derivs_t, jakob_t, bjakob_t, steper_state_t
  use Burn_data, ONLY : bn_meshMe
#include "constants.h"

  implicit none

!! argument declarations
  type(steper_state_t), intent(IN OUT) :: state
  integer, intent(IN) :: n
  real, intent(IN)    :: dydx(n), yscal(n), ratdum(:), htry, eps, btemp
  real, intent(INOUT) :: x, y(n)
  real, intent(OUT)   :: hdid, hnext
  procedure(derivs_t) :: derivs
  procedure(jakob_t) :: jakob
  procedure(bjakob_t) :: bjakob


!! local variables-----------------------------------------------
  integer, save       :: i,jtry
  integer, parameter  :: nmax=30, maxtry=400
  real, save          :: errmax,h,xsav,dysav(nmax),err(nmax),g1(nmax),  & 
       &                 g2(nmax),g3(nmax),g4(nmax),ysav(nmax),xx 

  real, parameter     ::  safety=0.9e0, grow=1.5e0, pgrow=-0.25e0,  & 
       &                  shrnk=0.5e0,  pshrnk=-1.0e0/3.0e0,errcon=0.1296e0

  real :: dydx_out(n)

!!  shampine parameter set 
  real, parameter     ::  gam =  1.0e0/2.0e0,    a21 =  2.0e0,  & 
       &                  a31 =  48.0e0/25.0e0,  a32 =  6.0e0/25.0e0,  & 
       &                  c21 = -8.0e0,          c31 =  372.0e0/25.0e0,  & 
       &                  c32 =  12.0e0/5.0e0,   c41 = -112.0e0/125.0e0,  & 
       &                  c42 = -54.0e0/125.0e0, c43 = -2.0e0/5.0e0,  & 
       &                  b1  =  19.0e0/9.0e0,   b2  =  1.0e0/2.0e0,  & 
       &                  b3  =  25.0e0/108.0e0, b4  =  125.0e0/108.0e0,  & 
       &                  e1  =  17.0e0/54.0e0,  e2  =  7.0e0/36.0e0,  & 
       &                  e3  =  0.0e0,          e4  =  125.0e0/108.0e0,  & 
       &                  c1x =  1.0e0/2.0e0,    c2x = -3.0e0/2.0e0,  & 
       &                  c3x =  121.0e0/50.0e0, c4x =  29.0e0/250.0e0,  & 
       &                  a2x =  1.0e0,          a3x = 3.0e0/5.0e0


!!  for the ma28 package
  integer, parameter  :: naij=200, n5 = 5*nmax, n8=8*nmax
  integer, parameter  :: nDummy = 1
  logical, save       :: firstCall = .false.
  integer, save       :: iloc(naij),jloc(naij), flag, nzo,  &
       &                 ivect(naij),jvect(naij),ikeep(n5),iw(n8)
  real, save          :: dfdy(naij,nDummy),amat(naij),w(nmax),u

!!----------------------------------------------------------------------


!!  get and copy the nonzero locations 
  if (.NOT. firstCall) then 
     firstCall = .true.

     nzo     = 0
     do i=1,naij
        iloc(i) = 0
        jloc(i) = 0
     end do

   !!  get the sparse pattern
     call bjakob(iloc,jloc,nzo,naij)

   !!  copy the location
     do i=1,nzo
        ivect(i) = iloc(i)
        jvect(i) = jloc(i)
     enddo

   !!  force the diagonal to be the pivot elements 
     do i=1,nzo 
        amat(i) = 1.0e-10 
        if (ivect(i) .eq. jvect(i)) amat(i) = 1.0e0 
     enddo

     u  = 0.1e0
     call ma28ad(n,nzo,amat,naij,iloc,naij,jloc,u,ikeep,iw,w,flag)

     if (flag .lt. 0) then
        if (bn_meshMe .EQ. MASTER_PE) print *, 'ERROR negative flag in ma28ad flag',flag
        if (bn_meshMe .EQ. MASTER_PE) print *, 'more than stpmax steps required in bn_roseMa28'
        call Logfile_stampMessage( '[bn_rosenMa28] negative flag returned from ma28ad')
        call Logfile_stampMessage( '[bn_rosenMa28] more than stpmax steps required')
        call Driver_abort('ERROR: more than stpmax steps required in bn_rosenMa28')
     end if

  end if  ! of firstCall =.false.


!!  store the initial values 
  xsav = x 
  do i=1,n 
     ysav(i)  = y(i) 
     dysav(i) = dydx(i) 
  enddo

!!  get the sparse jacobian in sparse_dfdy 
  call jakob(xsav,ysav,btemp,ratdum,dfdy,nzo,nDummy) 

!!  main loop 
  h = htry 
  do jtry = 1,maxtry 

   !!   form the a matrix and decompose it 
     xx = 1.0e0/(gam*h)
     do i=1,nzo 
        amat(i) = -dfdy(i,nDummy) 
        if (ivect(i) .eq. jvect(i)) amat(i) = xx + amat(i) 
     enddo

   !!  numeric decomp 
     call ma28bd(n,nzo,amat,naij,ivect,jvect,jloc,ikeep,iw,w,flag) 

     if (flag .lt. 0) then 
        if (bn_meshMe .EQ. MASTER_PE) print *, 'ERROR negative flag in ma28bd flag',flag
        call Logfile_stampMessage( '[bn_rosenMa28] negative flag returned from ma28bd')
        call Driver_abort('ERROR: negative return flag in ma28bd')
     end if


   !!  set up and solve the right hand side for g1 
     do i=1,n 
        g1(i) = dysav(i) 
     enddo
     call ma28cd(n,amat,naij,jloc,ikeep,g1,w,1) 


   !!  compute intermediate values of y,x and dydx 
     do i=1,n 
        y(i) = ysav(i) + a21 * g1(i) 
     enddo
     x = xsav + a2x * h 
     call derivs(x,y,btemp,ratdum,dydx_out) 


   !!  set up and solve the right hand side for g2 
     do i=1,n 
        g2(i) = dydx_out(i) + c21*g1(i)/h 
     enddo
     call ma28cd(n,amat,naij,jloc,ikeep,g2,w,1) 


   !!  compute intermediate values of y,x and dydx_out 
     do i=1,n 
        y(i) = ysav(i) + a31*g1(i) + a32*g2(i) 
     enddo
     x = xsav + a3x*h 
     call derivs(x,y,btemp,ratdum,dydx_out) 

   !!  set up and solve the right hand side for g3 
     do i=1,n 
        g3(i)  = dydx_out(i) + (c31*g1(i) + c32*g2(i))/h 
     enddo
     call ma28cd(n,amat,naij,jloc,ikeep,g3,w,1) 


   !!  set up and solve the right hand side for g4 
     do i=1,n 
        g4(i)  = dydx_out(i) + (c41*g1(i) + c42*g2(i) + c43*g3(i))/h 
     end do
     call ma28cd(n,amat,naij,jloc,ikeep,g4,w,1) 

   !!  compute the third and fourth order estimates of y 
     do i=1,n 
        y(i)   = ysav(i) + b1*g1(i) + b2*g2(i) + b3*g3(i) + b4*g4(i) 
        err(i) = e1*g1(i) + e2*g2(i) + e3*g3(i) + e4*g4(i) 
     enddo
     x = xsav + h 

     if (x .eq. xsav) then 
        if (bn_meshMe .EQ. MASTER_PE) print *, 'step size not significant in bn_rosenMa28'
        call Logfile_stampMessage( '[bn_rosenMa28] Step size not significant!')
        call Driver_abort('ERROR: step size not significant in bn_rosenMa28')
     end if


   !!  determine the scaled accuracy 
     errmax = 0.0e0 
     do i=1,n 
        errmax = max(errmax,abs(err(i)/yscal(i))) 
     enddo
     errmax = errmax/eps 

   !!  if the step succeded, compute the size of the next step and return 
     if (errmax .le. 1.0) then 
        hdid = h 
        if (errmax .gt. errcon) then 
           hnext = safety * h * errmax**pgrow 
        else 
           hnext = grow * h 
        end if
        return 

      !!   if the step did not succeed cut the stepsize and try again 
     else 
        hnext = safety * h * errmax**pshrnk 
        h     = sign(max(abs(hnext),shrnk*abs(h)),h) 
     end if
  enddo

!!  too many tries
  if (bn_meshMe .EQ. MASTER_PE) print *, 'ERROR exceeded maxtry in routine bn_rosenMa28' 
  call Logfile_stampMessage('[bn_rosenMa28] ERROR exceeded maxtry')
  call Driver_abort('ERROR: exceeded maxtry in bn_rosenMa28')

end subroutine bn_rosenMa28


