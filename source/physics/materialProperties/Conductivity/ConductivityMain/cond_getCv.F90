!!****if* source/physics/materialProperties/Conductivity/ConductivityMain/cond_getCv
!!
!! NAME
!!
!!  cond_getCv
!!
!! SYNOPSIS
!!
!!  cv = cond_getCv(real   (IN),OPTIONAL  :: xden,
!!                  real   (IN),OPTIONAL  :: xtemp,
!!                  integer(IN),OPTIONAL  :: component,
!!                  real   (IN),OPTIONAL,dimension(:)  :: massFrac,
!!                  real   (IN),OPTIONAL,dimension(:)  :: solnVec)
!!
!! DESCRIPTION
!!
!!   Get the appropriate value of specific heat capacity at constant
!!   volume, cv, for use by various physics units.
!!
!!   Currently used and implemented only inside the Conductivity unit,
!!   but the implementation has no dependencies on internals of the
!!   Conductivity unit and should be applicable more generally.
!!
!! ARGUMENTS
!!
!!   xden :  mass density, usually from DENS_VAR
!!
!!   xtemp : temperature
!!
!!   component : can select ion (1) or electron (2) component
!!
!!   massFrac : mass fractions
!!
!!   solnVec : full solution state vector for one cell
!!
!! NOTES
!!
!!  Technically, all arguments are optional. However, for proper functioning
!!  one of the following cases must apply:
!!    (a) at least the solnVec argument is present and has valid data;
!!    (b) at least the xden and xtemp arguments and addtionally, if required
!!        by the Eos implementation, the massFrac argument are present and
!!        have valid data.
!!
!!  If (a) applies, then xden, xtemp, and massFrac are ignored.
!!
!!  If possible, calls should be made as for case (a). That is, solnVec and,
!!  if applicable, component should be the only arguement provided in the
!!  function invocation, if the calling routine has a full state array
!!  available to pass as solnVec.
!!
!!***

real function cond_getCv(xden, xtemp, component, massFrac, solnVec)
  use Eos_interface, ONLY: Eos, Eos_getAbarZbar, Eos_getTempData
  implicit none

#include "Simulation.h"
#include "constants.h"
#include "Eos.h"

  real,   intent(IN),OPTIONAL :: xden
  real,   intent(IN),OPTIONAL :: xtemp
  integer,intent(IN),OPTIONAL :: component
  real   ,intent(IN),OPTIONAL :: massFrac(:)
  real   ,intent(IN),OPTIONAL :: solnVec(:)

  logical :: haveFullState
  integer :: whichCv, eCvIon, eCvEle
  real    :: abar,zbar
  real    :: cvOut
  real, dimension(EOS_NUM) :: eos_arr

  logical, dimension(EOS_VARS+1:EOS_NUM) :: mask

#ifdef TELE_VAR
  integer,parameter :: eosMode = MODE_DENS_TEMP_GATHER
#else
  integer,parameter :: eosMode = MODE_DENS_TEMP
#endif

  mask = .false.
  mask(EOS_CV)  = .true.
  mask(EOS_DET) = .true.
  whichCv = EOS_CV
#ifdef EOS_CVION
  eCvIon = EOS_CVION
#else
  eCvIon = EOS_CV
#endif
#ifdef EOS_CVELE
  eCvEle = EOS_CVELE
#else
  eCvEle = EOS_CV
#endif

  if (present(component)) then
     select case (component)
     case(1)
        whichCv = eCvIon
     case(2)
        whichCv = eCvEle
     end select
  end if
  mask(whichCv)  = .true.


  haveFullState = present(solnVec)

  if (haveFullState) then
     call Eos_getAbarZbar(solnVec=solnVec,abar=abar,zbar=zbar)
     eos_arr(EOS_DENS) = solnVec(DENS_VAR)
     call Eos_getTempData(solnVec,eos_arr,eosMode)
     eos_arr(EOS_ABAR) = abar
     eos_arr(EOS_ZBAR) = zbar

     if (NSPECIES > 0) then
        call Eos(eosMode,1,eos_arr,solnVec(SPECIES_BEGIN:SPECIES_END),mask)
     else
        call Eos(eosMode,1,eos_arr,mask=mask)
     end if

  else
     call Eos_getAbarZbar(abar=abar,zbar=zbar,massFrac=massFrac)
     eos_arr(EOS_DENS) = xden
     eos_arr(EOS_TEMP) = xtemp
#ifdef TELE_VAR
#   ifdef EOS_TEMPION
        eos_arr(EOS_TEMPION) = xtemp
#   endif
#   ifdef EOS_TEMPELE
        eos_arr(EOS_TEMPELE) = xtemp
#   endif
#   ifdef EOS_TEMPRAD
        eos_arr(EOS_TEMPRAD) = xtemp
#   endif
#endif

     eos_arr(EOS_ABAR) = abar
     eos_arr(EOS_ZBAR) = zbar

     if (NSPECIES > 0) then
        call Eos(eosMode,1,eos_arr,massFrac,mask)
     else
        call Eos(eosMode,1,eos_arr,mask=mask)
     end if

  end if

  zbar = eos_arr(EOS_ZBAR)
  cvOut = eos_arr(whichCv)
#ifdef TELE_VAR
  if (present(component)) then
     select case (component)
     case(1)
        if (whichCv == EOS_CV) then
           cvOut = cvOut / (1.0+zbar)
        end if
     case(2)
        if (whichCv == EOS_CV) then
           if (zbar .NE. 0.0) &
                cvOut = cvOut * zbar / (1.0+zbar)
        end if
     end select
  end if
#endif

  cond_getCv = cvOut

end function cond_getCv
