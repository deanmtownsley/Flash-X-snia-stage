!!****if* source/physics/sourceTerms/Burn/BurnIntegrate/bn_baderGift.F90
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
!!   bn_baderGift
!!
!! SYNOPSIS
!!  subroutine bn_baderMa28(real(IN)       ::y(:),
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
!! NOTES
!!  this file contains 2 routines to work with bader-duflhard integration steps
!!  for the integration nuclear reaction networks with 2 choices for the linear algebra.
!!
!!  routine bn_baderGift drives a bader-deuflhard step with gift algebra
!!  routine bn_baderStepGift is a bader-deuflhard stepper with ma28 algebra
!!
!!
!!***
!!---------------------------------------------------------------------------------

subroutine bn_baderGift(state,y,dydx,ratdum,nv,x,btemp,htry,eps,yscal,hdid,hnext, &
     &                       derivs,jakob,bjakob)

  use Burn_data, ONLY: aion
  use Driver_interface, ONLY : Driver_abort
  !  Ummm.... bit of a mystery why I can use the interfaces in bn_netIntRosen but not here.
  use bnIntegrate_interface, ONLY: bn_baderStepGift, bn_pzExtr, pzExtr_state_t
  use bnNetwork_interface, ONLY: derivs_t, jakob_t, bjakob_t, steper_state_t

  implicit none

!!  for dense analytic jacobians, gift algebra
!!
!!  semi-implicit extrapolation step for stiff ode's with monitoring
!!  of local truncation error to adjust stepsize. inputs are the dependent
!!  variable vector y(nv) and its derivative dydx(nv) at the starting of the
!!  independent variable x. also input are the stepsize to be attempted htry,
!!  the required accuracy eps, and the vector yscal against which the error is
!!  scaled. on output, y and x are replaced by their new values, hdid is the
!!  stepsize actually accomplished, and hnext is the estimated next stepsize.
!!  dervs is a user supplied function that computes the right hand side of
!!  the equations.
!!
!!  1/scalmx is the maximum increase in the step size allowed


!!  declare arguments
  type(steper_state_t), intent(IN OUT) :: state
  integer, intent(IN) :: nv
  real, intent(IN)    :: dydx(nv), yscal(nv), ratdum(:), htry, eps, btemp
  real, intent(INOUT) :: x, y(nv)
  real, intent(OUT)   :: hdid, hnext
  procedure(derivs_t) :: derivs
  procedure(jakob_t) :: jakob
  procedure(bjakob_t) :: bjakob

  integer, parameter :: nmax = 30, kmaxx=7, imax=kmaxx+1

  integer :: i,ii,iq,k,kk,km
  real :: errmax,fact,h,red,scale,xwork,wrkmin, &
          xest,xnew,err(kmaxx),yerr(nmax), &
          ysav(nmax),yseq(nmax), &
          dfdy(nmax,nmax)
  real, parameter :: safe1 = 0.25e0, safe2 = 0.7e0, redmax=1.0e-5, &
                     redmin = 0.7e0, tiny = 1.0e-30, scalmx = 0.1e0

  type(pzExtr_state_t) :: pzExtr_state

  logical :: converged

  associate(first => state%first, &
            reduct => state%reduct, &
            nvold => state%nvold, &
            eps1 => state%eps1, &
            epsold => state%epsold, &
            kmax => state%kmax, &
            kopt => state%kopt, &
            nseq => state%nseq, &
            a => state%a, &
            alf => state%alf)

!!  a new tolerance or a new number , so reinitialize
  if (eps .ne. epsold  .or.  nv .ne. nvold) then
     hnext = -1.0e29
     xnew  = -1.0e29
     eps1  = safe1 * eps

   !!  compute the work coefficients a_k
     a(1)  = nseq(1) + 1
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo

   !!  compute alf(k,q)
     do iq=2,kmaxx
        do k=1,iq-1
           alf(k,iq) = eps1**((a(k+1) - a(iq+1)) / &
                &               ((a(iq+1) - a(1) + 1.0e0) * (2*k + 1)))
        enddo
     enddo
     epsold = eps
     nvold  = nv

   !!  add cost of jacobians to work coefficients
     a(1)   = nv + a(1)
     do k=1,kmaxx
        a(k+1) = a(k) + nseq(k+1)
     enddo

   !!  determine optimal row number for convergence
     do kopt=2,kmaxx-1
        if (a(kopt+1) .gt. a(kopt)*alf(kopt-1,kopt)) exit
     enddo
     kmax = kopt
  end if


!!  save the starting values
  h    = htry
  do i=1,nv
     ysav(i)  = y(i)
  enddo

!!  get the dense jacobian in dfdy
  !      call Timers_start ("jacobian (burn)")

  call jakob(x,y,btemp,ratdum,dfdy,nv,nmax)

  !      call Timers_stop ("jacobian (burn)")

!!  a new stepsize or a new integration, re-establish the order window
  if (h .ne. hnext  .or.  x .ne. xnew) then
     first = .true.
     kopt = kmax
  end if
  reduct = .false.

  converged = .false.

outer_loop: do
  do i=1,nv
     y(i) = ysav(i)
  enddo

  ! Initialize red to avoid using uninitialized variable
  red = 1.0

  do k=1,kmax
     xnew = x + h

     if (xnew .eq. x) then
        write(*,*) 'stepsize too small in routine baderGift'
        write(*,110) xnew,x,k
110     format(1x,2e11.3,i6)

        do ii=1,nv
           write(*,111) y(ii),y(ii)*aion(ii),aion(ii),ii
111        format(1x,3e11.3,i6)
        enddo

        call Driver_abort('ERROR in baderGift: stepsize too small')
     end if

   !!  the semi implicit midpoint rules for this sequence
     call bn_baderStepGift(btemp,ysav,dydx,dfdy,ratdum,nmax,nv,x,h,nseq(k),yseq,derivs)

   !!  extrapolate the error to zero
     xest = (h/nseq(k))**2
     call bn_pzExtr(pzExtr_state,k,xest,yseq,y,yerr,nv)


   !!  compute normalized error estimate
     if (k .ne. 1) then
        errmax = tiny
        do i=1,nv
           errmax = max(errmax,abs(yerr(i)/yscal(i)))
        enddo
        errmax   = errmax/eps
        km = k - 1
        err(km) = (errmax/safe1)**(1.0e0/(2*km+1))
     end if

   !!  in order window
     if (k .ne. 1  .and. (k .ge. kopt-1  .or. first)) then

      !!   converged
        if (errmax .lt. 1.0) then
           converged = .true.
           exit  ! Exit the inner 'do k' loop
        end if


      !!   possible step size reductions
        if (k .eq. kmax  .or.  k .eq. kopt + 1) then
           red = safe2/err(km)
           red    = min(red,redmin)
           red    = max(red,redmax)
           h      = h * red
           reduct = .true.
           cycle outer_loop
        else if (k .eq. kopt) then
           if (alf(kopt-1,kopt) .lt. err(km)) then
              red = 1.0e0/err(km)
              red    = min(red,redmin)
              red    = max(red,redmax)
              h      = h * red
              reduct = .true.
              cycle outer_loop
           end if
        else if (kopt .eq. kmax) then
           if (alf(km,kmax-1) .lt. err(km)) then
              red = alf(km,kmax-1) * safe2/err(km)
              red    = min(red,redmin)
              red    = max(red,redmax)
              h      = h * red
              reduct = .true.
              cycle outer_loop
           end if
        else if (alf(km,kopt) .lt. err(km)) then
           red = alf(km,kopt-1)/err(km)
           red    = min(red,redmin)
           red    = max(red,redmax)
           h      = h * red
           reduct = .true.
           cycle outer_loop
        end if
     end if

  end do  ! End of 'do k=1,kmax' loop

  if (.not. converged) then
     red    = redmin
     h      = h * red
     reduct = .true.
     cycle outer_loop
  else
     exit outer_loop
  end if

end do outer_loop

!!  successful step; get optimal row for convergence and corresponding stepsize
  x = xnew
  hdid = h
  first = .false.
  wrkmin = 1.0e35
  do kk=1,km
     fact = max(err(kk),scalmx)
     xwork = fact * a(kk+1)
     if (xwork .lt. wrkmin) then
        scale  = fact
        wrkmin = xwork
        kopt   = kk + 1
     end if
  enddo

!!  check for possible order increase, but not if stepsize was just reduced
  hnext = h/scale
  if (kopt .ge. k  .and.  kopt .ne. kmax  .and.  .not.reduct) then
     fact = max(scale/alf(kopt-1,kopt),scalmx)
     if (a(kopt+1)*fact .le. wrkmin) then
        hnext = h/fact
        kopt = kopt + 1
     end if
  end if

  end associate

  return
end subroutine bn_baderGift

!!---------------------------------------------------------------------------------


subroutine bn_baderStepGift(btemp,y,dydx,dfdy,ratdum,nmax,n,xs,htot,nnstep,yout,  &
     &                      derivs)

  use bn_interface, ONLY: bn_gift
  use bnNetwork_interface, ONLY: derivs_t

  implicit none

!!
!!  an implicit midpoint stepper, for ma28 sparse linear algebra.
!!
!!  declare arguments
  integer, intent(IN) :: n, nmax
  integer, intent(IN) :: nnstep
  real, intent(IN)    :: btemp, xs, htot, y(n), dydx(n), dfdy(nmax,nmax), ratdum(:)
  real, intent(INOUT) :: yout(n)
  procedure(derivs_t) :: derivs

!!  declare local variables
  integer, parameter  :: nmaxx=30
  integer       ::  i,j,nn
  real           ::  h,x,del(nmaxx),ytemp(nmaxx)

!!  for the gift linear algebra
  integer, parameter  :: nmaxp1 = nmaxx + 1
  real          :: dmat(nmaxx,nmaxx), av(nmaxx,nmaxp1)


!!  stepsize this trip, and make the a matrix
  h = htot/nnstep

  !      call Timers_start ("algebra (burn)")

  do j=1,n
     do i=1,n
        dmat(i,j) = -h * dfdy(i,j)
     enddo
  enddo
  do i=1,n
     dmat(i,i) = 1.0e0 + dmat(i,i)
  end do

  !      call Timers_stop ("algebra (burn)")

!!  use yout as temporary storage; the first step
  do i=1,n
     yout(i) = h * dydx(i)
  enddo


  !      call Timers_start ("algebra (burn)")

  do j = 1, n
     do i = 1, n
        av(i,j) = dmat(i,j)
     enddo
  enddo
  do i = 1,n
     av(i,n+1) = yout(i)
  enddo
  call bn_gift(av,nmaxx,nmaxp1)
  do i = 1, n
     yout(i) = av(i,n+1)
  enddo

  !      call Timers_stop ("algebra (burn)")

  do i=1,n
     del(i)   = yout(i)
     ytemp(i) = y(i) + del(i)
  enddo
  x = xs + h

  !      call Timers_start ("derivs (burn)")

  call derivs(x,ytemp,btemp,ratdum,yout)

  !      call Timers_stop ("derivs (burn)")


!!  use yout as temporary storage; general step
  do nn=2,nnstep
     do i=1,n
        yout(i) = h*yout(i) - del(i)
     enddo

     !       call Timers_start ("algebra (burn)")

     do j = 1, n
        do i = 1, n
           av(i,j) = dmat(i,j)
        enddo
     enddo
     do i = 1,n
        av(i,n+1) = yout(i)
     enddo
     call bn_gift(av,nmaxx,nmaxp1)
     do i = 1, n
        yout(i) = av(i,n+1)
     enddo

     !       call Timers_stop ("algebra (burn)")


     do i=1,n
        del(i)   = del(i) + 2.0e0 * yout(i)
        ytemp(i) = ytemp(i) + del(i)
     enddo
     x = x + h


     !       call Timers_start ("derivs (burn)")

     call derivs(x,ytemp,btemp,ratdum,yout)

     !       call Timers_stop ("derivs (burn)")

  enddo

!!  take the last step
  do i=1,n
     yout(i) = h * yout(i) - del(i)
  enddo


  !      call Timers_start ("algebra (burn)")

  do j = 1, n
     do i = 1, n
        av(i,j) = dmat(i,j)
     enddo
  enddo
  do i = 1,n
     av(i,n+1) = yout(i)
  enddo
  call bn_gift(av,nmaxx,nmaxp1)
  do i = 1, n
     yout(i) = av(i,n+1)
  enddo

  !      call Timers_stop ("algebra (burn)")

  do i=1,n
     yout(i) = ytemp(i) + yout(i)
  enddo

  return

end subroutine bn_baderStepGift
