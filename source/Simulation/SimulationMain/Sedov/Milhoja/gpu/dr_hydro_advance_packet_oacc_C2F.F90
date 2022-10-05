!> @copyright Copyright 2022 UChicago Argonne, LLC and contributors
!!
!! @licenseblock
!! Licensed under the Apache License, Version 2.0 (the "License");
!! you may not use this file except in compliance with the License.
!!
!! Unless required by applicable law or agreed to in writing, software
!! distributed under the License is distributed on an "AS IS" BASIS,
!! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
!! See the License for the specific language governing permissions and
!! limitations under the License.
!! @endlicenseblock
!!
!! @file

#include "Milhoja.h"

#ifndef MILHOJA_OPENACC_OFFLOADING
#error "This file should only be compiled if using OpenACC offloading"
#endif

!> @details
!! A function that reinterprets the given C variables plucked from a data packet
!! into Flash-X compatible Fortran variables that are then passed onto the
!! automatically-generated Fortran routine that runs the actual work on the
!! data items in the associated packet.  This reinterpretation hinges on a
!! struct-of-arrays layout of the data packet.
!!
!! @todo The runtime should confirm at initialization that the casts will
!! always be correct for the given runtime configuration.
!! @todo Are the shapes of the different data blocks (e.g., U, face[XYZ]
!! minimized for the application)?
subroutine dr_hydro_advance_packet_oacc_c2f(         &
                    C_packet_h,                      &
                    C_dataQ_h,                       &
                    C_queue2_h,                      &
                    C_queue3_h,                      &
                    C_nTiles_h,                      &
                    C_nxbGC_h, C_nybGC_h, C_nzbGC_h, &
                    C_nCcVar_h, C_nFluxVar_h,        &
                    C_nTiles_d, C_dt_d,              &
                    C_deltas_start_d,                &
                    C_lo_start_d,   C_hi_start_d,    &
                    C_loGC_start_d,                  &
                    C_U_start_d,                     &
                    C_auxC_start_d,                  &
                    C_faceX_start_d,                 &
                    C_faceY_start_d,                 &
                    C_faceZ_start_d) bind(c)
    use iso_c_binding,     ONLY : C_PTR, &
                                  C_F_POINTER
    use openacc,           ONLY : acc_handle_kind
    use milhoja_types_mod, ONLY : MILHOJA_INT

    use dr_hydroAdvance_bundle_mod, ONLY : dr_hydroAdvance_packet_gpu_oacc

    implicit none

    type(C_PTR),          intent(IN), value :: C_packet_h
    integer(MILHOJA_INT), intent(IN), value :: C_dataQ_h
    integer(MILHOJA_INT), intent(IN), value :: C_queue2_h
    integer(MILHOJA_INT), intent(IN), value :: C_queue3_h
    integer(MILHOJA_INT), intent(IN), value :: C_nTiles_h
    integer(MILHOJA_INT), intent(IN), value :: C_nxbGC_h
    integer(MILHOJA_INT), intent(IN), value :: C_nybGC_h
    integer(MILHOJA_INT), intent(IN), value :: C_nzbGC_h
    integer(MILHOJA_INT), intent(IN), value :: C_nCcVar_h
    integer(MILHOJA_INT), intent(IN), value :: C_nFluxVar_h
    type(C_PTR),          intent(IN), value :: C_nTiles_d
    type(C_PTR),          intent(IN), value :: C_dt_d
    type(C_PTR),          intent(IN), value :: C_deltas_start_d
    type(C_PTR),          intent(IN), value :: C_lo_start_d
    type(C_PTR),          intent(IN), value :: C_hi_start_d
    type(C_PTR),          intent(IN), value :: C_loGC_start_d
    type(C_PTR),          intent(IN), value :: C_U_start_d
    type(C_PTR),          intent(IN), value :: C_auxC_start_d
    type(C_PTR),          intent(IN), value :: C_faceX_start_d
    type(C_PTR),          intent(IN), value :: C_faceY_start_d
    type(C_PTR),          intent(IN), value :: C_faceZ_start_d

    integer(kind=acc_handle_kind) :: F_dataQ_h
    integer(kind=acc_handle_kind) :: F_queue2_h
    integer(kind=acc_handle_kind) :: F_queue3_h
    integer                       :: F_nTiles_h
    integer                       :: F_nxbGC_h, F_nybGC_h, F_nzbGC_h
    integer                       :: F_nCcVar_h
    integer                       :: F_nFluxVar_h

    integer, pointer :: F_nTiles_d
    real,    pointer :: F_dt_d
    real,    pointer :: F_deltas_all_d(:,:)
    integer, pointer :: F_lo_all_d(:,:)
    integer, pointer :: F_hi_all_d(:,:)
    integer, pointer :: F_loGC_all_d(:,:)
    real,    pointer :: F_U_all_d(:,:,:,:,:)
    real,    pointer :: F_auxC_all_d(:,:,:,:)
    real,    pointer :: F_faceX_all_d(:,:,:,:,:)
    real,    pointer :: F_faceY_all_d(:,:,:,:,:)
    real,    pointer :: F_faceZ_all_d(:,:,:,:,:)

    !!!!!----- HOST-SIDE DATA TYPE MANAGEMENT
    ! Explicitly cast to Fortran default kinds without checking correctness
    ! We are assuming that the runtime has confirmed already that
    !   * MILHOJA_INT  = integer default kind
    !   * MILHOJA_REAL = real    default kind.
    F_dataQ_h    = INT(C_dataQ_h,  kind=acc_handle_kind)
    F_queue2_h   = INT(C_queue2_h, kind=acc_handle_kind)
    F_queue3_h   = INT(C_queue3_h, kind=acc_handle_kind)
    F_nTiles_h   = INT(C_nTiles_h)
    F_nxbGC_h    = INT(C_nxbGC_h)
    F_nybGC_h    = INT(C_nybGC_h)
    F_nzbGC_h    = INT(C_nzbGC_h)
    F_nCcVar_h   = INT(C_nCcVar_h)
    F_nFluxVar_h = INT(C_nFluxVar_h)

!    write(*,*) "[C2F] OpenAcc Queue = ", F_dataQ_h
!    write(*,*) "[C2F] nTiles        = ", F_nTiles_h
!    write(*,*) "[C2F] nxbGC         = ", F_nxbGC_h
!    write(*,*) "[C2F] nybGC         = ", F_nybGC_h
!    write(*,*) "[C2F] nzbGC         = ", F_nzbGC_h
!    write(*,*) "[C2F] nCcVar        = ", F_nCcVar_h
!    write(*,*) "[C2F] nFluxVar      = ", F_nFluxVar_h

    !!!!!----- DEVICE-SIDE DATA TYPE MANAGEMENT
    ! Get Fortran view of C device pointers
    ! Reshaping dictated by internal struct-of-arrays layout of data packet
    CALL C_F_POINTER(C_nTiles_d,       F_nTiles_d)
    CALL C_F_POINTER(C_dt_d,           F_dt_d)
    CALL C_F_POINTER(C_deltas_start_d, F_deltas_all_d, shape=[3, F_nTiles_h])
    CALL C_F_POINTER(C_lo_start_d,     F_lo_all_d,     shape=[3, F_nTiles_h])
    CALL C_F_POINTER(C_hi_start_d,     F_hi_all_d,     shape=[3, F_nTiles_h])
    CALL C_F_POINTER(C_loGC_start_d,   F_loGC_all_d,   shape=[3, F_nTiles_h])
    CALL C_F_POINTER(C_U_start_d,      F_U_all_d,      shape=[F_nxbGC_h,   F_nybGC_h,   F_nzbGC_h,   F_nCcVar_h,   F_nTiles_h])
    CALL C_F_POINTER(C_auxC_start_d,   F_auxC_all_d,   shape=[F_nxbGC_h,   F_nybGC_h,   F_nzbGC_h,                 F_nTiles_h])
    CALL C_F_POINTER(C_faceX_start_d,  F_faceX_all_d,  shape=[F_nxbGC_h+1, F_nybGC_h,   F_nzbGC_h,   F_nFluxVar_h, F_nTiles_h])
#if MILHOJA_NDIM == 1
    CALL C_F_POINTER(C_faceY_start_d,  F_faceY_all_d,  shape=[1,           1,           1,           1,            F_nTiles_h])
#else
    CALL C_F_POINTER(C_faceY_start_d,  F_faceY_all_d,  shape=[F_nxbGC_h,   F_nybGC_h+1, F_nzbGC_h,   F_nFluxVar_h, F_nTiles_h])
#endif
#if MILHOJA_NDIM <= 2
    CALL C_F_POINTER(C_faceZ_start_d,  F_faceZ_all_d,  shape=[1,           1,           1,           1,            F_nTiles_h])
#else
    CALL C_F_POINTER(C_faceZ_start_d,  F_faceZ_all_d,  shape=[F_nxbGC_h,   F_nybGC_h,   F_nzbGC_h+1, F_nFluxVar_h, F_nTiles_h])
#endif

    !!!!!----- PASS DATA TO TARGET-SPECIFIC STATIC FORTRAN LAYER
    CALL dr_hydroAdvance_packet_gpu_oacc(C_packet_h,                        &
                                         F_dataQ_h, F_queue2_h, F_queue3_h, &
                                         F_nTiles_d, F_dt_d,                &
                                         F_deltas_all_d,                    &
                                         F_lo_all_d,   F_hi_all_d,          &
                                         F_loGC_all_d,                      &
                                         F_U_all_d,                         &
                                         F_auxC_all_d,                      &
                                         F_faceX_all_d,                     &
                                         F_faceY_all_d,                     &
                                         F_faceZ_all_d)
end subroutine dr_hydro_advance_packet_oacc_c2f

