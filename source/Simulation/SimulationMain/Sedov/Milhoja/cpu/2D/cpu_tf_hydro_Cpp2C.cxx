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
#include <Milhoja_IntVect.h>
#include <Milhoja_RealVect.h>
#include <Milhoja_Tile.h>
#include <Milhoja_interface_error_codes.h>

#include "Tile_cpu_tf_hydro.h"

extern "C" {
    //----- C DECLARATION OF FORTRAN ROUTINE WITH C-COMPATIBLE INTERFACE
    void cpu_tf_hydro_c2f (
        const milhoja::Real external_hydro_op1_dt,
        const int external_hydro_op1_eosMode,
        const void* tile_deltas,
        const void* tile_hi,
        const void* tile_interior,
        const void* tile_lo,
        const void* CC_1,
        const void* scratch_hydro_op1_auxC,
        const void* scratch_hydro_op1_flX,
        const void* scratch_hydro_op1_flY,
        const void* scratch_hydro_op1_flZ,
        const void* lbdd_CC_1,
        const void* lbdd_scratch_hydro_op1_auxC,
        const void* lbdd_scratch_hydro_op1_flX,
        const void* lbdd_scratch_hydro_op1_flY
    );

    //----- C DECLARATION OF ACTUAL TASK FUNCTION TO PASS TO RUNTIME
    void  cpu_tf_hydro_cpp2c(const int threadIndex,
                             milhoja::DataItem* dataItem) {
        Tile_cpu_tf_hydro*  wrapper = dynamic_cast<Tile_cpu_tf_hydro*>(dataItem);
        milhoja::Tile*      tileDesc = wrapper->tile_.get();
    
        const milhoja::IntVect   tile_lo     = tileDesc->lo();
        const milhoja::IntVect   tile_hi     = tileDesc->hi();
        const milhoja::IntVect   tile_lbound = tileDesc->loGC();
        const milhoja::RealVect  tile_deltas = tileDesc->deltas();
        milhoja::Real*           CC_1        = tileDesc->dataPtr();

        milhoja::Real* scratch_hydro_op1_auxC = 
                 static_cast<milhoja::Real*>(Tile_cpu_tf_hydro::scratch_hydro_op1_auxc_)
                + Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_AUXC_SIZE_ * threadIndex;
   
        milhoja::Real* scratch_hydro_op1_flX = 
                 static_cast<milhoja::Real*>(Tile_cpu_tf_hydro::scratch_hydro_op1_flx_)
                + Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_FLX_SIZE_ * threadIndex;
        milhoja::Real* scratch_hydro_op1_flY = 
                 static_cast<milhoja::Real*>(Tile_cpu_tf_hydro::scratch_hydro_op1_fly_)
                + Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_FLY_SIZE_ * threadIndex;
        milhoja::Real* scratch_hydro_op1_flZ = 
                 static_cast<milhoja::Real*>(Tile_cpu_tf_hydro::scratch_hydro_op1_flz_)
                + Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_FLZ_SIZE_ * threadIndex;

        milhoja::Real  tile_deltas_array[]
                = {tile_deltas.I(), tile_deltas.J(), tile_deltas.K()};
        int   tile_lo_array[] = {tile_lo.I(), tile_lo.J(), tile_lo.K()};
        int   tile_hi_array[] = {tile_hi.I(), tile_hi.J(), tile_hi.K()};
        // Create as 1D with ordering setup for mapping on 2D Fortran array
        int   tile_interior[] = {tile_lo.I(), tile_hi.I(),
                                 tile_lo.J(), tile_hi.J(),
                                 tile_lo.K(), tile_hi.K()};
        int   lbdd_scratch_hydro_op1_auxC[] = {tile_lo.I() - MILHOJA_K1D,
                                               tile_lo.J() - MILHOJA_K2D,
                                               tile_lo.K() - MILHOJA_K3D};
        // No need for variable masking in tile wrappers since the full data
        // array is in the memory system and available for use.  Always start
        // with first variable?  If we have to use variable masking (See C2F
        // layer), then we would need to set CC_1 to the start variable here.
        int   lbdd_CC_1[]
               = {tile_lbound.I(), tile_lbound.J(), tile_lbound.K(), 1};
        int   lbdd_scratch_hydro_op1_flX[]
               = {tile_lo.I(), tile_lo.J(), tile_lo.K(), 1};
        int   lbdd_scratch_hydro_op1_flY[]
               = {tile_lo.I(), tile_lo.J(), tile_lo.K(), 1};

        cpu_tf_hydro_c2f(wrapper->external_hydro_op1_dt_,
                         wrapper->external_hydro_op1_eosMode_,
                         static_cast<void*>(tile_deltas_array),
                         static_cast<void*>(tile_hi_array),
                         static_cast<void*>(tile_interior),
                         static_cast<void*>(tile_lo_array),
                         static_cast<void*>(CC_1),
                         static_cast<void*>(scratch_hydro_op1_auxC),
                         static_cast<void*>(scratch_hydro_op1_flX),
                         static_cast<void*>(scratch_hydro_op1_flY),
                         static_cast<void*>(scratch_hydro_op1_flZ),
                         static_cast<void*>(lbdd_CC_1),
                         static_cast<void*>(lbdd_scratch_hydro_op1_auxC),
                         static_cast<void*>(lbdd_scratch_hydro_op1_flX),
                         static_cast<void*>(lbdd_scratch_hydro_op1_flY));
    }
}
