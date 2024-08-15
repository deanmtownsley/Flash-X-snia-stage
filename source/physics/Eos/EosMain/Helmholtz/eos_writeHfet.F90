!> @copyright Copyright 2023 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!!   Licensed under the Apache License, Version 2.0 (the "License");
!!   you may not use this file except in compliance with the License.
!!
!!   Unless required by applicable law or agreed to in writing, software
!!   distributed under the License is distributed on an "AS IS" BASIS,
!!   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!!   See the License for the specific language governing permissions and
!!   limitations under the License.
!! @endlicenseblock
!!
!! @file
!> @ingroup physics_Eos
!!
!! @brief  write coefficients for tabular helmholtz in a binary file
!!
!! @param  n -- number of variables in the arrays
!! @param  f --  Helmholtz free energy
!!


subroutine eos_writeHfet(n, f)
   use Driver_interface, ONLY : Driver_abort

   implicit none

#include "Simulation.h"

   integer, intent(in) :: n
   real, intent(in) :: f(n)
!! Local variables
!!  fileUnit is a file variable
   integer, parameter ::  fileUnit = 36
   integer          :: numRead, ioStat



#ifdef DEBUG
  if (n<0) then
    call Driver_abort("[eos_writeHfet]  n must be positive")
  endif 
#endif

!! Open the file
  open (fileUnit,FILE='helm_table.bdat',ACTION='WRITE',STATUS='REPLACE',FORM='UNFORMATTED',IOSTAT=ioStat)
    if (ioStat .NE. 0)     call Driver_abort("[eos_writeHfet]  file open failure!")

!! Start writing

  write(fileUnit,ERR=101) f
  
!! close up and return
  close (fileUnit,IOSTAT=ioStat)
   if (ioStat .NE. 0)  call Driver_abort("[eos_writeHfet]  couldn't close file!")
  
  return

!!  Abort statements

101    call Driver_abort("[eos_writeHfet]  failed write on f!")

end subroutine eos_writeHfet


