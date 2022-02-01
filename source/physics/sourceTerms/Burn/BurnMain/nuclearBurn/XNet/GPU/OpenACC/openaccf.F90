module openaccf
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
   use, intrinsic :: iso_c_binding
   implicit none

   interface
      subroutine acc_map_data(hostptr,devptr,bytes) &
            bind(c,name="acc_map_data")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: hostptr
         type(c_ptr), value :: devptr
         integer(c_size_t), value :: bytes
      end subroutine acc_map_data

      subroutine acc_unmap_data(hostptr) &
            bind(c,name="acc_unmap_data")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: hostptr
      end subroutine acc_unmap_data

      type(c_ptr) function acc_deviceptr(hostptr) &
            bind(c,name="acc_deviceptr")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: hostptr
      end function acc_deviceptr

      type(c_ptr) function acc_hostptr(devptr) &
            bind(c,name="acc_hostptr")
         use, intrinsic :: iso_c_binding
         type(c_ptr), value :: devptr
      end function acc_hostptr
   end interface

end module