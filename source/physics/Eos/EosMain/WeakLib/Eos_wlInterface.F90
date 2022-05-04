Module Eos_nucInterface

  implicit none

  interface
    subroutine Eos_wlOneZone(xDens,xTemp,xYe,xEner,xPres,xEntr,xdedt,xCs2,xXp,xXn,xXa,xXh,xAbar,xVar,varID,mode)
    implicit none
    real, intent(INOUT) :: xDens, xYe
    real, intent(INOUT) :: xTemp, xEner, xEntr, xPres
    integer, intent(IN) :: mode, varID
    real, intent(OUT) :: xXp, xXn, xXa,xXh,xdedt,xCs2,xVar,xAbar
    real :: xZbar,xMu_e,xMu_n,xMu_p,xMuhat
    end subroutine Eos_wlOneZone
  end interface

  interface
     subroutine Eos_wlDetectBounce (postBounce,bounceTime,centralDens,centralEntr)
       implicit none
       logical, intent(OUT) :: postBounce
       real, optional, intent(OUT) :: bounceTime, centralDens, centralEntr
     end subroutine Eos_wlDetectBounce
  end interface

  interface
     subroutine Eos_wlEnerShift (energyShift)
       implicit none
       real, intent(OUT) :: energyShift
     end subroutine Eos_wlEnerShift
  end interface

end Module Eos_nucInterface
  
