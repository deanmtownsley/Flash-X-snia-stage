!   Returns the properties of the unburned fuel and first-stage ashes based on
!   the initial abundances as passed in
!
! subroutine argements:
!
!   xc12initial    :  initial abundance of C12 in this cell
!   xne20initial   :  initial abundance of Ne20 in this cell
!   xne22initial   :  initial abundance of Ne22 in this cell
!   ye_f, ye_a     :  electron fraction for fuel and ash (electrons/Baryon)
!   yi_f, yi_a     :  Ion "mole fraction" for fuel and ash (ions/Baryon)
!   qbar_f, qbar_a :  Average binding energy per baryon for fuel and ash
!
!  Dean Townsley 2008
!

subroutine bn_paraFuelAshProperties(xc12initial, xne20initial, xne22initial, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

  use Burn_data, ONLY :  ye_c12,  yi_c12,  q_c12, &
                           ye_o16,  yi_o16,  q_o16, &
                           ye_ne22, yi_ne22, q_ne22,&
                           ye_mg24, yi_mg24, q_mg24,&
                           ye_ne20, yi_ne20, q_ne20

  implicit none
  real, intent(in)  :: xc12initial, xne20initial, xne22initial
  real, intent(out) :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a
  real :: xo16initial

! does it still hold that the ash is just burning to mg24 if ne20 is fuel? Well, the ne20 and c12 burn together so this is more or less okay.
! Alternately we could say the ash is an admixture of mg24 and the primary products of ne20 burning depending on xne20cone.

  ! simple implementation for now
  !  ash is just result of C12 & Ne20->Mg24
  !  likely to eventually be based on some manner of table interpolation
  xo16initial = 1.0-xc12initial-xne20initial-xne22initial

  ye_f = xc12initial*ye_c12 + xo16initial*ye_o16 + xne20initial*ye_ne20 + xne22initial*ye_ne22
  ye_a = (xc12initial+xne20initial)*ye_mg24 + xne22initial*ye_ne22 + xo16initial*ye_o16

  yi_f =  xc12initial*yi_c12 + xo16initial*yi_o16 + xne20initial*yi_ne20 + xne22initial*yi_ne22
  yi_a = (xc12initial+xne20initial)*yi_mg24 + xne22initial*yi_ne22 + xo16initial*yi_o16

  qbar_f =  xc12initial*q_c12 + xo16initial*q_o16 + xne20initial*q_ne20 + xne22initial*q_ne22
  qbar_a = (xc12initial+xne20initial)*q_mg24 + xne22initial*q_ne22 + xo16initial*q_o16

end subroutine
