!!****ih* source/physics/materialProperties/Conductivity/localAPI/cond_interface
!! NAME
!!  Conductivity_interface
!!
!! SYNOPSIS
!!  USE cond_interface, ONLY: cond_getCv
!!
!! DESCRIPTION
!!  This is an internal interface module for the Conductivity Unit that defines some
!!  private interfaces.
!!***

Module cond_interface

interface
   real function cond_getCv(xden, xtemp, component, massFrac, solnVec)
     implicit none
     real,   intent(IN),OPTIONAL :: xden
     real,   intent(IN),OPTIONAL :: xtemp
     integer,intent(IN),OPTIONAL :: component
     real   ,intent(IN),OPTIONAL :: massFrac(:)
     real   ,intent(IN),OPTIONAL :: solnVec(:)
   end function cond_getCv
end interface

end Module cond_interface
