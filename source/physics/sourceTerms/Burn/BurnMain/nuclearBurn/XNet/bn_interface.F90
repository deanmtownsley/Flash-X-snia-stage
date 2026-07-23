!!****ih* source/physics/sourceTerms/Burn/BurnMain/nuclearBurn/XNet/bn_interface
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
!! SYNOPSIS
!!   use bn_interface
!!
!! DESCRIPTION
!!
!! This is the header file for the Burn module that defines its
!! private interfaces.
!!
!!***

Module bn_interface

#include "Simulation.h"
#include "constants.h"

  interface
     subroutine bn_xnetInit(data_dir,data_desc)
        implicit none
        character(*), intent(in) :: data_dir
        character(80), intent(out) :: data_desc
     end subroutine bn_xnetInit
  end interface

  interface
     subroutine bn_xnetFinalize()
        implicit none
     end subroutine bn_xnetFinalize
  end interface

  interface
     subroutine bn_burner(tstep,temp,density,xIn,xOut,sdotRate,burnedZone,zone,kstep)
       implicit none
       logical, intent(IN), dimension(:)              :: burnedZone
       integer, intent(IN), dimension(:,:)            :: zone
       real, intent(IN)                               :: tstep
       real, intent(IN), dimension(:)                 :: temp,density
       real, intent(OUT), dimension(size(temp))       :: sdotRate
       real, intent(IN), dimension(:,:)               :: xIn
       real, intent(OUT), dimension(size(xIn,1),size(xIn,2)) :: xOut
       integer, intent(OUT)                           :: kstep
     end subroutine bn_burner
  end interface

  interface
     subroutine bn_azbar(xmass, ymass, abar, zbar, z2bar, ytot1, bye)
       implicit none
       real, intent(IN) :: xmass(NSPECIES)
       real, intent(OUT) :: ymass(NSPECIES)
       real, intent(OUT) :: abar, zbar, z2bar, ytot1, bye
     end subroutine bn_azbar
  end interface

  interface
     subroutine bn_ecapnuc(etakep,temp,rpen,rnep,spen,snep)
       implicit none
       real, intent(IN)     :: etakep, temp
       real, intent(OUT)    :: rpen,rnep,spen,snep
     end subroutine bn_ecapnuc
  end interface

  interface
     real function bn_ifermi12(f)
       implicit none
       real, intent(IN) :: f
     end function bn_ifermi12
  end interface

  interface
     subroutine bn_mazurek(btemp,bden,y56,ye,rn56ec,sn56ec)
       implicit none
       real, intent(IN) :: btemp, bden
       real, intent(IN) :: y56, ye
       real, intent(OUT):: rn56ec,sn56ec

     end subroutine bn_mazurek
  end interface

  interface
     subroutine bn_mcord(i,j,iloc,jloc,nzo,np,eloc,nterm,np2)
       implicit none
       integer, intent(IN)  ::  i,j, np, np2
       integer, intent(INOUT) :: nterm, nzo
       integer, intent(OUT) ::  iloc(np),jloc(np),eloc(np2)
     end subroutine bn_mcord
  end interface

  type :: screen4_state_t
     real :: qlam0z
     real :: gamp
     real :: taufac
     real :: gamef
     real :: tau12
     real :: alph12
     real :: h12w
     real :: h12
     real :: xlgfac
     real :: cc
     real :: xx
     real :: gamp14
     real :: alp123
     real :: xni
     real :: aa
     real :: bb
     real :: dd
     real :: btempi
     real :: btemp_old
     real :: den_old
     real :: zbarr_old
     real :: abarr_old
  end type screen4_state_t

  interface
     subroutine bn_screen4(state, zbarr, abarr, z2barr, z1, a1, z2, a2, &
                           jscreen, init, &
                           btemp, bden, zs13, zs13inv, zhat, zhat2, lzav, aznut, &
                           scorr)
       import :: screen4_state_t
       implicit none
       type(screen4_state_t), intent(INOUT) :: state
       integer, intent(IN) :: jscreen
       logical, intent(IN) :: init
       real, intent(IN) :: abarr, zbarr, z2barr, z1, a1, z2, a2
       real, intent(IN) :: btemp, bden
       real, intent(OUT) :: scorr
       real, intent(INOUT) :: zs13(:), zs13inv(:), zhat(:), zhat2(:), lzav(:), aznut(:)   ! size of nrat
     end subroutine bn_screen4
  end interface

  interface
     subroutine bn_sneutx(btemp,bden,abar,zbar, &
                          sneut)
       implicit none
       real, intent(IN) :: btemp,bden,abar,zbar
       real, intent(OUT) :: sneut
     end subroutine bn_sneutx
  end interface

end Module bn_interface
