!!
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
!!
!!**
subroutine Stencils_advectWeno2d(rhs,phi,u,v,dx,dy,ix1,ix2,jy1,jy2,center,facex,facey)

  implicit none
  !----Imported variables
  real, dimension(:,:,:), intent(inout):: rhs
  real, dimension(:,:,:), intent(in) :: phi,u,v
  real, intent(in) :: dx,dy
  integer, intent(in) :: ix1,ix2,jy1,jy2
  integer, intent(in) :: center,facex,facey
        
  !------Local variables
  real ::  err,eps,d0,d1,d2,eyl,eyr,ur,ul,vr,vl,& 
           s1r,s2r,s3r,s4r,s5r,s1l,s2l,s3l,s4l,s5l, &
           rIS1r,rIS2r,rIS3r,rIS1l,rIS2l,rIS3l, &
           aT1r,aT2r,aT3r,aT1l,aT2l,aT3l, &
           a1r,a2r,a3r,a1l,a2l,a3l, &
           fT1r,fT2r,fT3r,fT1l,fT2l,fT3l, &
           frx,flx,fry,fly
  integer :: i,j,k
 
  k = 1
  eps = 1e-15

  do j=jy1,jy2
     do i=ix1,ix2

     ul = center*u(i,j,k) + &
          facex*0.5*(u(i,j,k) + u(i-1,j,k)) + &
          facey*0.5*(u(i,j,k) + u(i,j-1,k))

     ur = center*u(i+1,j,k) + &
          facex*0.5*(u(i+1,j,k)+u(i,j,k)) + &
          facey*0.5*(u(i+1,j,k)+u(i+1,j-1,k))

     vl = center*v(i,j,k) + &
          facex*0.5*(v(i,j,k) + v(i-1,j,k)) + &
          facey*0.5*(v(i,j,k) + v(i,j-1,k))

     vr = center*v(i,j+1,k) + &
          facex*0.5*(v(i,j+1,k) + v(i-1,j+1,k)) + &
          facey*0.5*(v(i,j+1,k) + v(i,j,k))

     !______________________Advection Terms_______________________!

     !----------------- WENO3 X-Direction ------------!
     if (ur .gt. 0) then     ! u = (+) Downwind

        s1r = phi(i-2,j,k)
        s2r = phi(i-1,j,k)
        s3r = phi(i,j,k)
        s4r = phi(i+1,j,k)
        s5r = phi(i+2,j,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 1./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 3./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
        fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
        fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

      else                  ! u = (-) Upwind

        s1r = phi(i-1,j,k)
        s2r = phi(i,j,k)
        s3r = phi(i+1,j,k)
        s4r = phi(i+2,j,k)
        s5r = phi(i+3,j,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 3./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 1./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
        fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
        fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

     end if

     if (ul .gt. 0) then     ! u = (+) Downwind  

        s1l = phi(i-3,j,k)
        s2l = phi(i-2,j,k)
        s3l = phi(i-1,j,k)
        s4l = phi(i,j,k)
        s5l = phi(i+1,j,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 1./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 3./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
        fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
        fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

     else                   ! u = (-) Upwind

        s1l = phi(i-2,j,k)
        s2l = phi(i-1,j,k)
        s3l = phi(i,j,k)
        s4l = phi(i+1,j,k)
        s5l = phi(i+2,j,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 3./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 1./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
        fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
        fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

     end if

     !---------------------------------------------------------
     !- WENO3 interpolated PHI values at cell face
     !---------------------------------------------------------
     frx = a1r*fT1r + a2r*fT2r + a3r*fT3r
     flx = a1l*fT1l + a2l*fT2l + a3l*fT3l
     !---------------------------------------------------------
     !---------------------------------------------------------

     !----------------- WENO3 Y-Direction ------------!
     if (vr .gt. 0) then     ! u = (+) Downwind

        s1r = phi(i,j-2,k)
        s2r = phi(i,j-1,k)
        s3r = phi(i,j,k)
        s4r = phi(i,j+1,k)
        s5r = phi(i,j+2,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 1./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 3./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r =  2./6.*s1r - 7./6.*s2r + 11./6.*s3r
        fT2r = -1./6.*s2r + 5./6.*s3r +  2./6.*s4r
        fT3r =  2./6.*s3r + 5./6.*s4r -  1./6.*s5r

     else                   ! u = (-) Upwind

        s1r = phi(i,j-1,k)
        s2r = phi(i,j,k)
        s3r = phi(i,j+1,k)
        s4r = phi(i,j+2,k)
        s5r = phi(i,j+3,k)

        rIS1r = 13./12.*(    s1r  - 2.*s2r +    s3r )**2. &
              +  1./4. *(    s1r  - 4.*s2r + 3.*s3r )**2.
        rIS2r = 13./12.*(    s2r  - 2.*s3r +    s4r )**2. &
              +  1./4. *(    s2r           -    s4r )**2.
        rIS3r = 13./12.*(    s3r  - 2.*s4r +    s5r )**2. &
              +  1./4. *( 3.*s3r  - 4.*s4r +    s5r )**2.

        aT1r = 3./10. /  ( eps + rIS1r )**2.
        aT2r = 6./10. /  ( eps + rIS2r )**2.
        aT3r = 1./10. /  ( eps + rIS3r )**2.

        a1r = aT1r / ( aT1r + aT2r +aT3r )
        a2r = aT2r / ( aT1r + aT2r +aT3r )
        a3r = aT3r / ( aT1r + aT2r +aT3r )

        fT1r = -1./6.*s1r + 5./6.*s2r +  2./6.*s3r
        fT2r =  2./6.*s2r + 5./6.*s3r -  1./6.*s4r
        fT3r =  11./6.*s3r - 7./6.*s4r + 2./6.*s5r

     end if

     if (vl .gt. 0) then     ! u = (+) Downwind

        s1l = phi(i,j-3,k)
        s2l = phi(i,j-2,k)
        s3l = phi(i,j-1,k)
        s4l = phi(i,j,k)
        s5l = phi(i,j+1,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 1./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 3./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l =  2./6.*s1l - 7./6.*s2l + 11./6.*s3l
        fT2l = -1./6.*s2l + 5./6.*s3l +  2./6.*s4l
        fT3l =  2./6.*s3l + 5./6.*s4l -  1./6.*s5l

     else                    ! u = (-) Upwind

        s1l = phi(i,j-2,k)
        s2l = phi(i,j-1,k)
        s3l = phi(i,j,k)
        s4l = phi(i,j+1,k)
        s5l = phi(i,j+2,k)

        rIS1l = 13./12.*(    s1l  - 2.*s2l +    s3l )**2. &
              +  1./4. *(    s1l  - 4.*s2l + 3.*s3l )**2.
        rIS2l = 13./12.*(    s2l  - 2.*s3l +    s4l )**2. &
              +  1./4. *(    s2l           -    s4l )**2.
        rIS3l = 13./12.*(    s3l  - 2.*s4l +    s5l )**2. &
              +  1./4. *( 3.*s3l  - 4.*s4l +    s5l )**2.

        aT1l = 3./10. /  ( eps + rIS1l )**2.
        aT2l = 6./10. /  ( eps + rIS2l )**2.
        aT3l = 1./10. /  ( eps + rIS3l )**2.

        a1l = aT1l / ( aT1l + aT2l +aT3l )
        a2l = aT2l / ( aT1l + aT2l +aT3l )
        a3l = aT3l / ( aT1l + aT2l +aT3l )

        fT1l = -1./6.*s1l + 5./6.*s2l +  2./6.*s3l
        fT2l =  2./6.*s2l + 5./6.*s3l -  1./6.*s4l
        fT3l =  11./6.*s3l - 7./6.*s4l + 2./6.*s5l

     end if

     !---------------------------------------------------------
     !- WENO3 interpolated PHI values at cell face
     !---------------------------------------------------------
     fry = a1r*fT1r + a2r*fT2r + a3r*fT3r
     fly = a1l*fT1l + a2l*fT2l + a3l*fT3l
     !---------------------------------------------------------
     !---------------------------------------------------------

!_______________________________RHS TERM______________________________________!

     rhs(i,j,k) = rhs(i,j,k) - (frx*ur  - flx*ul)/dx  - (fry*vr  - fly*vl)/dy &
                             + phi(i,j,k)*(ur - ul)/dx + phi(i,j,k)*(vr - vl)/dy

    end do
  end do 

  return

end subroutine Stencils_advectWeno2d
