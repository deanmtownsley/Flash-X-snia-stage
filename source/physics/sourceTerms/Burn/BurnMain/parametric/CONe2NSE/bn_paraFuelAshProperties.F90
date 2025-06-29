!   Returns the properties of the unburned fuel and first-stage ashes based on
!   the initial abundances as passed in
!
! subroutine argements:
!
!   xc12initial    :  initial abundance of C12 in this cell
!   xne22initial   :  initial abundance of Ne22 in this cell
!   ye_f, ye_a     :  electron fraction for fuel and ash (electrons/Baryon)
!   yi_f, yi_a     :  Ion "mole fraction" for fuel and ash (ions/Baryon)
!   qbar_f, qbar_a :  Average binding energy per baryon for fuel and ash
!
!  Dean Townsley 2008
!

subroutine bn_paraFuelAshProperties(xc12initial, xne22initial, ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a)

  use Burn_data, ONLY :  ye_c12,  yi_c12,  q_c12, &
                           ye_o16,  yi_o16,  q_o16, &
                           ye_ne22, yi_ne22, q_ne22,&
                           ye_mg24, yi_mg24, q_mg24

  implicit none
  real, intent(in)  :: xc12initial, xne22initial
  real, intent(out) :: ye_f, ye_a, yi_f, yi_a, qbar_f, qbar_a

  ! simple implementation for now
  !  ash is just result of C12->Mg24
  !  likely to eventually be based on some manner of table interpolation
  ye_f = xc12initial*ye_c12 + xne22initial*ye_ne22 + (1.0-xc12initial-xne22initial)*ye_o16
  ye_a = xc12initial*ye_mg24 + xne22initial*ye_ne22 + (1.0-xc12initial-xne22initial)*ye_o16

  yi_f = xc12initial*yi_c12 + xne22initial*yi_ne22 + (1.0-xc12initial-xne22initial)*yi_o16
  yi_a = xc12initial*yi_mg24 + xne22initial*yi_ne22 + (1.0-xc12initial-xne22initial)*yi_o16

  qbar_f = xc12initial*q_c12 + xne22initial*q_ne22 + (1.0-xc12initial-xne22initial)*q_o16
  qbar_a = xc12initial*q_mg24 + xne22initial*q_ne22 + (1.0-xc12initial-xne22initial)*q_o16

end subroutine
