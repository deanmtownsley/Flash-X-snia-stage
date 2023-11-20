/**
 * @copyright Copyright 2022 UChicago Argonne, LLC and contributors
 *
 * @licenseblock
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * @endlicenseblock
 *
 * @file
 */

#include <Milhoja.h>
#include <Milhoja_real.h>
#include <Milhoja_interface_error_codes.h>

#include "Tile_cpu_tf_hydro.h"

extern "C" {
    //----- C DECLARATION OF FORTRAN ROUTINE WITH C-COMPATIBLE INTERFACE
    int instantiate_hydro_advance_wrapper_c(const milhoja::Real dt, void** wrapper) {
        if        ( wrapper == nullptr) {
            std::cerr << "[instantiate_hydro_advance_wrapper_c] wrapper is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        } else if (*wrapper != nullptr) {
            std::cerr << "[instantiate_hydro_advance_wrapper_c] *wrapper not NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_NOT_NULL;
        }

        try {
            *wrapper = static_cast<void*>(new Tile_cpu_tf_hydro{dt});
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_WRAPPER;
        } catch (...) {
            std::cerr << "[instantiate_hydro_advance_wrapper_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_WRAPPER;
        }

        return MILHOJA_SUCCESS;
    }

    int delete_hydro_advance_wrapper_c(void* wrapper) {
        if (wrapper == nullptr) {
            std::cerr << "[delete_hydro_advance_wrapper_c] wrapper is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        delete static_cast<Tile_cpu_tf_hydro*>(wrapper);

        return MILHOJA_SUCCESS;
    }

    int get_dt_wrapper_c(void* wrapper, milhoja::Real* dt) {
        if (!wrapper || !dt) {
            std::cerr << "[get_dt_wrapper_c] wrapper or dt NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }

        try {
            Tile_cpu_tf_hydro*   tile = static_cast<Tile_cpu_tf_hydro*>(wrapper);
            *dt = tile->dt_;
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_ACCESS_WRAPPER;
        } catch (...) {
            std::cerr << "[get_dt_wrapper_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_ACCESS_WRAPPER;
        }

        return MILHOJA_SUCCESS;
    }

    int acquire_scratch_wrapper_c(void) {
        Tile_cpu_tf_hydro::acquireScratch();

        return MILHOJA_SUCCESS;
    }

    int release_scratch_wrapper_c(void) {
        Tile_cpu_tf_hydro::releaseScratch();

        return MILHOJA_SUCCESS;
    }

    int get_scratch_auxc_wrapper_c(void* wrapper, const int threadIdx,
                                   milhoja::Real** auxC) {
        if (!wrapper || !auxC) {
            std::cerr << "[get_scratch_auxc_wrapper_c] wrapper or auxC NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        } else if (*auxC) {
            std::cerr << "[get_scratch_auxc_wrapper_c] *auxC not NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_NOT_NULL;
        } else if (threadIdx < 0) {
            std::cerr << "[get_scratch_auxc_wrapper_c] Negative threadIdx" << std::endl;
            return MILHOJA_ERROR_NEGATIVE_VALUE_FOR_UINT;
        }

        try {
            Tile_cpu_tf_hydro*   tile = static_cast<Tile_cpu_tf_hydro*>(wrapper);
            *auxC =   static_cast<milhoja::Real*>(tile->hydro_op1_auxc_)
                    + Tile_cpu_tf_hydro::HYDRO_OP1_AUXC_SIZE_ * threadIdx;
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_ACCESS_WRAPPER;
        } catch (...) {
            std::cerr << "[get_scratch_auxc_wrapper_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_ACCESS_WRAPPER;
        }

        return MILHOJA_SUCCESS;
    }
}
