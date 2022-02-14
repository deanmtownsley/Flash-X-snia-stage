subroutine ut_qsortInt(inOutArray, arrayLength, ascOrderArg)
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
  implicit none
  integer, dimension(:), intent(INOUT) :: inOutArray
  integer, intent(IN) :: arrayLength
  logical, optional, intent(IN) :: ascOrderArg
  logical :: ascOrder

  if (present(ascOrderArg)) then
     ascOrder = ascOrderArg
  else
     ascOrder = .true.
  end if

  if (ascOrder .eqv. .true.) then
     call ut_qsort_int_asc(inOutArray, arrayLength)
  else
     call ut_qsort_int_Desc(inOutArray, arrayLength)  
  end if

end subroutine ut_qsortInt


subroutine ut_qsortFloat(inOutArray, arrayLength, ascOrderArg)
  implicit none
  real*4, dimension(:), intent(INOUT) :: inOutArray
  integer, intent(IN) :: arrayLength
  logical, optional, intent(IN) :: ascOrderArg
  logical :: ascOrder

  if (present(ascOrderArg)) then
     ascOrder = ascOrderArg
  else
     ascOrder = .true.
  end if

  if (ascOrder .eqv. .true.) then
     call ut_qsort_float_asc(inOutArray, arrayLength)
  else
     call ut_qsort_float_desc(inOutArray, arrayLength)  
  end if

end subroutine ut_qsortFloat


subroutine ut_qsortDouble(inOutArray, arrayLength, ascOrderArg)
  implicit none
  real*8, dimension(:), intent(INOUT) :: inOutArray
  integer, intent(IN) :: arrayLength
  logical, optional, intent(IN) :: ascOrderArg
  logical :: ascOrder

  if (present(ascOrderArg)) then
     ascOrder = ascOrderArg
  else
     ascOrder = .true.
  end if

  if (ascOrder .eqv. .true.) then
     call ut_qsort_double_asc(inOutArray, arrayLength)
  else
     call ut_qsort_double_desc(inOutArray, arrayLength)  
  end if

end subroutine ut_qsortDouble
