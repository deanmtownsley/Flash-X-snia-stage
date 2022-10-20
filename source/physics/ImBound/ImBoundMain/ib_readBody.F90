!!***if* source/physics/ImBound/ImBoundMain/ib_readBody
!!
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
!!
!!***
subroutine ib_readBody(body, bodyFile)

   use ImBound_type, ONLY: ImBound_type_t
   use Driver_interface, ONLY: Driver_abort
   use HDF5

   implicit none
   type(ImBound_type_t), intent(inout) :: body
   character(len=*), intent(in) :: bodyFile

   !-------------------------------------------------------------------------------
   integer(HID_T)                 :: file
   integer                        :: h5err
   integer(HID_T)                 :: dset
   integer(HSIZE_T), dimension(3) :: dims

   call h5open_f(h5err)
   call h5fopen_f(trim(bodyFile), H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F
   if (h5err < 0) call Driver_abort('Unable to open body file')

   dims = (/1, 1, 1/)

   call h5dopen_f(file, "numElems", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read body/numElems')
   call h5dread_f(dset, H5T_NATIVE_INTEGER, body%numElems, dims, h5err)
   call h5dclose_f(dset, h5err)

   allocate (body%elems(body%numElems))

   dims = (/body%numElems, 1, 1/)

   call h5dopen_f(file, "elems/xA", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/xA')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%xA, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/yA", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/yA')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%yA, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/xB", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/xB')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%xB, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/yB", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/yB')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%yB, dims, h5err)
   call h5dclose_f(dset, h5err)

   call h5fclose_f(file, h5err)
   call h5close_f(h5err)
   
   body%dim = 2 !! need to read this from file 
end subroutine ib_readBody
