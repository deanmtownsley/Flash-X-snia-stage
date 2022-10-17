!!****if* source/physics/Hydro/HydroMain/Hydro_mapBcType
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
!!  Hydro_mapBcType
!!
!!
!!  For more details see the documentation of the NULL implementation
!!
!!***

#include "Simulation.h"

#ifdef DEBUG_ALL
#define DEBUG_GRID
#endif
#define DEBUG_GRID

subroutine Hydro_mapBcType(bcTypeToApply,bcTypeFromGrid,varIndex,gridDataStruct, &
     axis,face,idest)

#include "constants.h"
  
  use Driver_interface, ONLY : Driver_abort
  use Hydro_data, ONLY: hy_meshMe
  

  implicit none
  
  integer, intent(OUT) :: bcTypeToApply
  integer, intent(in) :: bcTypeFromGrid,varIndex,gridDataStruct,axis,face
  integer,intent(IN),OPTIONAL:: idest

  integer :: var,i,j,k,n,m,strt,fin, varCount,bcVecEnd
  logical :: validGridDataStruct
  logical :: isMagComponent


#ifdef DEBUG_GRID

  validGridDataStruct = .false.
  validGridDataStruct= (gridDataStruct == CENTER).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEX).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEY).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == FACEZ).or.validGridDataStruct
  validGridDataStruct= (gridDataStruct == WORK).or.validGridDataStruct
  
  if(.not.validGridDataStruct) then
     if (hy_meshMe .EQ. MASTER_PE) print *, "Hydro_mapBcType: gridDataStruct set to improper value"
     if (hy_meshMe .EQ. MASTER_PE) print *, "gridDataStruct must = CENTER,FACEX,FACEY,FACEZ,WORK " // &
          " (defined in constants.h)"
     call Driver_abort("Hydro_mapBcType gridDataStruct must be one of CENTER,FACEX,FACEY,FACEZ,WORK(see constants.h)")
  end if

  if((gridDataStruct==WORK).and.(varIndex/=1)) &
       call Driver_abort("Hydro_mapBcType: varCount be 1 for work array")
  if((gridDataStruct==CENTER).and.(varIndex > NUNK_VARS)) &
       call Driver_abort("Hydro_mapBcType: varIndex be <= NUNK_VARS for unk array")
       

#endif
  
  bcTypeToApply = 0
  return  
  isMagComponent = .FALSE.
#ifdef MAG_FACE_VAR
  if (varIndex==MAG_FACE_VAR) isMagComponent = .TRUE.
#endif
#ifdef MAGI_FACE_VAR
  if (varIndex==MAGI_FACE_VAR) isMagComponent = .TRUE.
#endif


  if (isMagComponent) then
                 !! Note: A special consideration is applied for face variables using outflow BC.
                 !! The conventional outflow treatment breaks divB=0 condition that is subject to
                 !! satisfy for the constrained-transport (CT) facevariables. 
                 !! There are two ways to do this depending on to which variables Neumann BC applies:
                 !! (1) apply Neumann condition to the two transverse (to the boundary or asix) fields,
                 !!     then calculate the normal field using the divB=0 relationship as a closure,
                 !! (2) apply Neumann condition to the normal field, then seek for the transverse fields
                 !!     that satisfy divB=0. A simple way is to set them as constants.
                 !!     For example, if we set Neumann for Bx along x-direction, then assigning 0 to
                 !!     By and Bz will satisfy divB=0 automatically.
                 !! We adopt the approach (2) for FLASH implementation as there are two advantages:
                 !! (i)  BC handlings in Paramesh can only work for each gridDataStruct at a time.
                 !!      This becomes a limitation for the approach (1) because we need
                 !!      FACEX, FACEY and FACEZ all at the same time to close Bx using Neumann By & Bz.
                 !!      On the other hand, the approach (2) only requires each FACE* gridDataStruct
                 !!      without needing the other two for closure.
                 !! (ii) The normal fields have more meaningful usages in calculations. For instance,
                 !!      They are directly used for n+1/2 Riemann states, whereas the corresponding
                 !!      states for transverse fields are always reconstructed. In this sense, it is
                 !!      better to apply Neumann for the normal fields, while any constant values can
                 !!      be assigned to the transverse fields. Such constant values are then filled,
                 !!      satisfying divB=0, without affecting Riemann state reconstructions because
                 !!      it is the (Neumann-treated) cell-centered B-fields that are used for 
                 !!      the transverse fields reconstructions.
     if (gridDataStruct == FACEX .OR. gridDataStruct == FACEY .OR. gridDataStruct == FACEZ) then
        if ((gridDataStruct == FACEX .and. axis == IAXIS) .or. &
            (gridDataStruct == FACEY .and. axis == JAXIS) .or. &
            (gridDataStruct == FACEZ .and. axis == KAXIS)) then
           select case (bcTypeFromGrid)
           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
              bcTypeToApply = REFLECTING
           end select
        end if
        if ((gridDataStruct == FACEX .and. axis .NE. IAXIS) .or. &
            (gridDataStruct == FACEY .and. axis .NE. JAXIS) .or. &
            (gridDataStruct == FACEZ .and. axis .NE. KAXIS)) then
           select case (bcTypeFromGrid)
           case(OUTFLOW,HYDROSTATIC_F2_NVOUT,HYDROSTATIC_F2_NVDIODE,HYDROSTATIC_F2_NVREFL, &
                HYDROSTATIC_NVOUT,HYDROSTATIC_NVDIODE,HYDROSTATIC_NVREFL)
              bcTypeToApply = GRIDBC_ZERO
           end select
        end if
     end if
  end if

  return
end subroutine Hydro_mapBcType








