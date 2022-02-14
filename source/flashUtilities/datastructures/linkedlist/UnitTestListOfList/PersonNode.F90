module PersonNode
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
  use Driver, ONLY : Driver_abort
  implicit none

  type person_node
     character (len=80) :: name
     integer :: age
     type(person_node), pointer :: next, prev
  end type person_node

contains

  !The user must provide: create_node, destroy_node and print_node.
  subroutine create_node(item)
    implicit none
    type(person_node), pointer :: item
    integer :: err    
    allocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abort ("[PersonNode::create_node]: "//&
            "Memory cannot be allocated")
    end if
    item % name = "NULL"
    item % age = -1
    nullify(item % next, item % prev)
  end subroutine create_node


  subroutine destroy_node(item)
    implicit none
    type(person_node), pointer  :: item
    integer :: err

    deallocate(item, STAT=err)
    if (err /= 0) then
       call Driver_abort ("[PersonNode::destroy_node]: "//&
            "Memory cannot be deallocated")
    end if
    nullify(item)
  end subroutine destroy_node


  subroutine print_node(unitNumber, item)    
    implicit none
    integer, intent(IN) :: unitNumber
    type(person_node), pointer :: item
    write(unitNumber,'(3A,i3)') "Name: ", trim(item % name), &
         ", age: ", item % age
  end subroutine print_node
end module PersonNode
