module UnitTest_NodeObject
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
  type node
     integer :: someInt, threadThatTouched
     type(node), pointer :: next, prev
  end type node

  contains

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(node), pointer :: item
    integer :: err    
    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abort &
            ("[create_node]: No more heap memory")
    end if
    item % SomeInt = -1
    item % threadThatTouched = -1
    nullify(item % next, item % prev)
  end subroutine create_node


  subroutine destroy_node(item)
    implicit none
    type(node), pointer  :: item    
    deallocate(item); nullify(item)
  end subroutine destroy_node


  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(node), pointer :: item
    write(unitNumber,*) "[print_node]: Value", item % SomeInt, &
         "(Optional Thread)", item % threadThatTouched
  end subroutine print_node


  !Just whack this here, so we can test without FLASH.
  subroutine Driver_abort(msg)
    implicit none
    character (len=*), intent(IN) :: msg
    print *, "ERROR!!!"
    print *, msg
    stop
  end subroutine Driver_abort

end module UnitTest_NodeObject
