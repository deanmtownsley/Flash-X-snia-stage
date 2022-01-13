## Lines starting with ## are comments inside template file
## All other lines including empty lines are non-comments
## 
## This file is a template for generating the setup_getFlashUnits.F90
## source file.  For syntax of this file see "Readme.template".
##
##
## VALID VARIABLE NAMES FOR THIS TEMPLATE
##
## unit_names -> list of strings containing names of units 
##
## Each unitname must have length atmost 80 chars, maybe.
##
!!****f* object/setup_flashUnits
!!
!! NAME
!!
!!  setup_getFlashUnits
!!
!!
!! SYNOPSIS
!!
!!
!!  call setup_getFlashUnits(character(len=MAX_STRING_LENGTH) :: unit_names(1:$num_units))
!!
!!  call setup_getFlashUnits(character(len=MAX_STRING_LENGTH) :: unit_names(1:NUM_UNITS))
!!
!!
!! DESCRIPTION
!!
!!  Return a character array of size NUM_UNITS containing
!!  the names of all of the FLASH units used to assemble
!!  the current executable
!!
!!  The unit_names variable should be declared as
!!
!!    include "constants.h"
!!  
!!    character (len=MAX_STRING_LENGTH) :: flash_units(numUnits)
!!
!!  or
!!
!!    character(len=MAX_STRING_LENGTH),ALLOCATABLE :: flash_units(:)
!!
!!
!!  The length of each character string is set to MAX_STRING_LENGTH,
!!  which is defined in the constants.h file.
!!  The proper number of elements to allocate, called numUnits above,
!!  can be inquired at run time with
!!
!!   call setup_getNumFlashUnits(numUnits)
!!
!!  so that an ALLOCATABLE array can be used for flash_units.
!!
!!***

  subroutine setup_getFlashUnits(unit_names)

#include "constants.h"
    implicit none

    integer, PARAMETER :: NUM_UNITS = %(COUNT_unit_names)s
    character (len=MAX_STRING_LENGTH) :: unit_names(NUM_UNITS)
    integer :: i

    i = 0

    i = i + 1; unit_names(i) = &
"%(unit_names!"\n    i = i + 1; unit_names(i) = &\n")s"


    return

  end subroutine setup_getFlashUnits

  subroutine setup_getNumFlashUnits(numUnits)

    implicit none

    integer, intent(out) :: numUnits
    integer, PARAMETER :: NUM_UNITS = %(COUNT_unit_names)s

    numUnits = NUM_UNITS

    return

  end subroutine setup_getNumFlashUnits

