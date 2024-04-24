!!****if* source/physics/Eos/EosMain/Helmholtz/eos_writeHfet
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
!!
!!  eos_writeHfet
!!
!! SYNOPSIS
!!
!!  eos_writeHfet(integer(IN) :: n, 
!!            real(OUT) :: f(n), 
!!            real(OUT) :: fd(n), 
!!            real(OUT) :: ft(n), 
!!            real(OUT) :: fdd(n), 
!!            real(OUT) :: ftt(n),
!!            real(OUT) :: fdt(n), 
!!            real(OUT) :: fddt(n), 
!!            real(OUT) :: fdtt(n), 
!!            real(OUT) :: fddtt(n), 
!!            real(OUT) :: dpdf(n), 
!!            real(OUT) :: dpdfd(n),
!!            real(OUT) :: dpdft(n), 
!!            real(OUT) :: dpdfdt(n), 
!!            real(OUT) :: ef(n), 
!!            real(OUT) :: efd(n), 
!!            real(OUT) :: eft(n), 
!!            real(OUT) :: efdt(n), 
!!            real(OUT) :: xf(n), 
!!            real(OUT) :: xfd(n), 
!!            real(OUT) :: xft(n), 
!!            real(OUT) :: xfdt(n)  )
!!
!! DESCRIPTION
!!
!!  Read binary data file containing coefficients for tabular helmholtz
!!
!! ARGUMENTS
!!
!!     n -- number of variables in the arrays
!!     f --  Helmholtz free energy
!!     fd --  derivative of f wrt density
!!     ft --  derivative of f wrt temperature
!!     fdd --  second derivative of f wrt density
!!     ftt --  second derivative of f wrt temperature
!!     fdt --  second derivative of f wrt density and temperature
!!     fddt --  third derivative of f wrt density^2 and temperature
!!     fdtt --  third derivative of f wrt density and temperature^2 e.g. dF/(dd)(dt^2)
!!     fddtt --  fourth derivative of f wrt density^2 and temperature^2
!!     dpdf --  pressure derivative 
!!     dpdfd -- 
!!     dpdft --  
!!     dpdfdt --  
!!     ef --  electron chemical potential
!!     efd --  
!!     eft --  
!!     efdt --  
!!     xf --  number density
!!     xfd --  
!!     xft --  
!!     xfdt --
!!
!!  NOTE
!! 
!!    See Timmes and Swesty, 2000, AJSS, "The Accuracy, Consistency, and Speed of an Electron-Positron
!!    Equation of State Based on Table Interpolation of the Helmholtz Free Energy"
!!
!!***

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


