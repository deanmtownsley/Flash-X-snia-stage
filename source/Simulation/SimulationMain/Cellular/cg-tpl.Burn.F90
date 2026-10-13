! This is a CG-Kit template file that can be inserted as a node
! in the ORCHA's recipe using FlashX_RecipeTools. This file can be
! inserted by using the `FlashX_RecipeTools.TimeStepRecipe.add_tpl` function
! to insert certain code lines in four different places of TimeAdvance.F90:
!    use_interface: before the "implicit none"
!                   usually for "use MODULE, only: SUBROUTINE" line
!    var_definition: after the "implicit none"
!                    usually for a variable declarations
!    var_initialization: after the variable declarations.
!                        usually for initializing variables
!    execute: main body of the TimeAdvance
!
! This example shows how to invoke a Burn subroutine *outside* of Milhoja.


!<_connector:use_interface>
use Burn_interface, ONLY: Burn


!<_connector:var_definition>
! None


!<_connector:var_initialization>
! None


!<_connector:execute>
call Timers_start("sourceTerms")
call Burn(dt)
call Timers_stop("sourceTerms")

