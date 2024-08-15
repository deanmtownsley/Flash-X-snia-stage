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
#include <Milhoja_TileWrapper.h>
#include <Milhoja_interface_error_codes.h>

extern "C" {
    //----- C DECLARATION OF FORTRAN ROUTINE WITH C-COMPATIBLE INTERFACE
    int instantiate_wrapper_c(void** wrapper) {
        if        ( wrapper == nullptr) {
            std::cerr << "[instantiate_wrapper_c] wrapper is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        } else if (*wrapper != nullptr) {
            std::cerr << "[instantiate_wrapper_c] *wrapper not NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_NOT_NULL;
        }

        try {
            *wrapper = static_cast<void*>(new milhoja::TileWrapper());
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_WRAPPER;
        } catch (...) {
            std::cerr << "[instantiate_wrapper_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_WRAPPER;
        }

        return MILHOJA_SUCCESS;
    }

    int delete_wrapper_c(void* wrapper) {
        if (wrapper == nullptr) {
            std::cerr << "[delete_wrapper_c] wrapper is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        delete static_cast<milhoja::TileWrapper*>(wrapper);

        return MILHOJA_SUCCESS;
    }
}

