module UnitTest_ListObject
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
  use UnitTest_NodeObject
  implicit none

  type list
     type(node), pointer :: H, T
  end type list

contains

  !An include file with the basic list operations AND 
  !a higher order function which operates on the nodes 
  !of the list.  The higher order function requires the 
  !following C preprocessor definition.
#define CPP_NODE_DEFINITION \
  use UnitTest_NodeObject
#include "ut_listMethods.includeF90"
#undef CPP_NODE_DEFINITION
end module UnitTest_ListObject
