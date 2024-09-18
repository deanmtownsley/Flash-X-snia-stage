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
   integer(HSIZE_T), dimension(3) :: dsetDims

   call h5open_f(h5err)
   call h5fopen_f(trim(bodyFile), H5F_ACC_RDONLY_F, file, h5err) !H5F_ACC_RDONLY_F
   if (h5err < 0) call Driver_abort('Unable to open body file')

   dsetDims = (/1, 1, 1/)

   call h5dopen_f(file, "numElems", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read body/numElems')
   call h5dread_f(dset, H5T_NATIVE_INTEGER, body%numElems, dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "dims", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read body/dims')
   call h5dread_f(dset, H5T_NATIVE_INTEGER, body%dims, dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   allocate (body%elems(body%numElems))

   dsetDims = (/body%numElems, 1, 1/)

   call h5dopen_f(file, "elems/xA", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/xA')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pA(1), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/yA", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/yA')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pA(2), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/xB", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/xB')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pB(1), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/yB", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/yB')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pB(2), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/xCenter", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/xCenter')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%center(1), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/yCenter", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/yCenter')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%center(2), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/xNorm", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/xNorm')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%normal(1), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   call h5dopen_f(file, "elems/yNorm", dset, h5err)
   if (h5err < 0) call Driver_abort('Unable to read elems/yNorm')
   call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%normal(2), dsetDims, h5err)
   call h5dclose_f(dset, h5err)

   body%boundBox(:, 1) = (/minval(body%elems(:)%center(1)), &
                           maxval(body%elems(:)%center(1))/)

   body%boundBox(:, 2) = (/minval(body%elems(:)%center(2)), &
                           maxval(body%elems(:)%center(2))/)

   if (body%dims == 3) then
      ! Add read statement for 3d body

      call h5dopen_f(file, "elems/zA", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/zA')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pA(3), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      call h5dopen_f(file, "elems/zB", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/zB')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pB(3), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      call h5dopen_f(file, "elems/zCenter", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/zCenter')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%center(3), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      call h5dopen_f(file, "elems/zNorm", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/zNorm')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%normal(3), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      call h5dopen_f(file, "elems/xC", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/xC')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pC(1), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      call h5dopen_f(file, "elems/yC", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/yC')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pC(2), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      call h5dopen_f(file, "elems/zC", dset, h5err)
      if (h5err < 0) call Driver_abort('Unable to read elems/zC')
      call h5dread_f(dset, H5T_NATIVE_DOUBLE, body%elems(:)%pC(3), dsetDims, h5err)
      call h5dclose_f(dset, h5err)

      body%boundBox(:, 3) = (/minval(body%elems(:)%center(3)), &
                              maxval(body%elems(:)%center(3))/)

   end if

   call h5fclose_f(file, h5err)
   call h5close_f(h5err)

end subroutine ib_readBody
