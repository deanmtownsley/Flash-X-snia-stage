

#include <iostream>
#include <Milhoja.h>
#include <Milhoja_real.h>
#include <Milhoja_interface_error_codes.h>

#include "DataPacket_gpu_tf_hydro.h"

#ifndef MILHOJA_OPENACC_OFFLOADING
#error "This file should only be compiled if using OpenACC offloading"
#endif

using milhoja::Real;

extern "C" {
    //----- C DECLARATION OF FORTRAN ROUTINE WITH C-COMPATIBLE INTERFACE
    void gpu_tf_hydro_c2f (
    void* packet_h,
    const int queue1_h,
    const int _nTiles_h,
    const void* _nTiles_d,
    const void* _external_hydro_op1_dt_d,
    const void* _tile_deltas_d,
    const void* _tile_lo_d,
    const void* _tile_hi_d,
    const void* _tile_lbound_d,
    const void* _CC_1_d,
    const void* _scratch_hydro_op1_auxC_d,
    const void* _scratch_hydro_op1_flX_d,
    const void* _scratch_hydro_op1_flY_d,
    const void* _scratch_hydro_op1_flZ_d
    
    );

    int instantiate_gpu_tf_hydro_packet_c (
    real external_hydro_op1_dt,void** packet
    
    ) {
        if ( packet == nullptr) {
            std::cerr << "[instantiate_gpu_tf_hydro_packet_c] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        } else if (*packet != nullptr) {
            std::cerr << "[instantiate_gpu_tf_hydro_packet_c] *packet not NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_NOT_NULL;
        }

        try {
            *packet = static_cast<void*>(new DataPacket_gpu_tf_hydro(
            external_hydro_op1_dt
            ));
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_PACKET;
        } catch (...) {
            std::cerr << "[instantiate_gpu_tf_hydro_packet_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_PACKET;
        }

        return MILHOJA_SUCCESS;
    }

    int delete_gpu_tf_hydro_packet_c (void* packet) {
        if (packet == nullptr) {
            std::cerr << "[delete_gpu_tf_hydro_packet_c] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        delete static_cast< DataPacket_gpu_tf_hydro *>(packet);
        return MILHOJA_SUCCESS;
    }

    int release_gpu_tf_hydro_extra_queue_c (void* packet, const int id) {
        std::cerr << "[release_gpu_tf_hydro_extra_queue_c] Packet does not have extra queues." << std::endl;
        return MILHOJA_ERROR_UNABLE_TO_RELEASE_STREAM;
    }

    //----- C TASK FUNCTION TO BE CALLED BY RUNTIME
    void gpu_tf_hydro_cpp2c (const int threadIndex, void* dataItem_h) {
        DataPacket_gpu_tf_hydro* packet_h = static_cast<DataPacket_gpu_tf_hydro*>(dataItem_h);
        const int queue1_h = packet_h->asynchronousQueue();
        const int _nTiles_h = packet_h->_nTiles_h;
        

        void* _nTiles_d = static_cast<void*>( packet_h->_nTiles_d );
        void* _external_hydro_op1_dt_d = static_cast<void*>( packet_h->_external_hydro_op1_dt_d );
        void* _tile_deltas_d = static_cast<void*>( packet_h->_tile_deltas_d );
        void* _tile_lo_d = static_cast<void*>( packet_h->_tile_lo_d );
        void* _tile_hi_d = static_cast<void*>( packet_h->_tile_hi_d );
        void* _tile_lbound_d = static_cast<void*>( packet_h->_tile_lbound_d );
        void* _CC_1_d = static_cast<void*>( packet_h->_CC_1_d );
        void* _scratch_hydro_op1_auxC_d = static_cast<void*>( packet_h->_scratch_hydro_op1_auxC_d );
        void* _scratch_hydro_op1_flX_d = static_cast<void*>( packet_h->_scratch_hydro_op1_flX_d );
        void* _scratch_hydro_op1_flY_d = static_cast<void*>( packet_h->_scratch_hydro_op1_flY_d );
        void* _scratch_hydro_op1_flZ_d = static_cast<void*>( packet_h->_scratch_hydro_op1_flZ_d );
        

        // Pass data packet info to C-to-Fortran Reinterpretation Layer
        gpu_tf_hydro_c2f (
        packet_h,
        queue1_h,
        _nTiles_h,
        _nTiles_d,
        _external_hydro_op1_dt_d,
        _tile_deltas_d,
        _tile_lo_d,
        _tile_hi_d,
        _tile_lbound_d,
        _CC_1_d,
        _scratch_hydro_op1_auxC_d,
        _scratch_hydro_op1_flX_d,
        _scratch_hydro_op1_flY_d,
        _scratch_hydro_op1_flZ_d
        
        );
    }
}