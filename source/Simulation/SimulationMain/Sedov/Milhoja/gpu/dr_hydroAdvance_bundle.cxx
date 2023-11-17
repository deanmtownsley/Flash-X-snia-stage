

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
    // -------------changed function name
    void dr_hydro_advance_packet_oacc_c2f (
    void* packet_h,
    const int queue1_h,
    const int queue2_h,
    const int queue3_h,
    const int _nTiles_h,
    const void* _nTiles_d,
    const void* _dt_d,
    const void* _tile_deltas_d,
    const void* _tile_lo_d,
    const void* _tile_hi_d,
    const void* _tile_lbound_d,
    const void* _U_d,
    const void* _hydro_op1_auxc_d,
    const void* _hydro_op1_flX_d,
    const void* _hydro_op1_flY_d,
    const void* _hydro_op1_flZ_d
    
    );

    // ---------------changed function name
    int instantiate_hydro_advance_packet_c (
    real dt,void** packet
    
        ) {
        if ( packet == nullptr) {
            std::cerr << "[instantiate_DataPacket_gpu_tf_hydro_C] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        } else if (*packet != nullptr) {
            std::cerr << "[instantiate_DataPacket_gpu_tf_hydro_C] *packet not NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_NOT_NULL;
        }

        try {
            *packet = static_cast<void*>(new DataPacket_gpu_tf_hydro(
            dt
            ));
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_PACKET;
        } catch (...) {
            std::cerr << "[instantiate_DataPacket_gpu_tf_hydro_C] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_PACKET;
        }

        return MILHOJA_SUCCESS;
    }

    // ----------------changed function name
    int delete_hydro_advance_packet_c (void* packet) {
        if (packet == nullptr) {
            std::cerr << "[delete_DataPacket_gpu_tf_hydro_C] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        delete static_cast< DataPacket_gpu_tf_hydro *>(packet);
        return MILHOJA_SUCCESS;
    }

    // ----------------changed function name
    int release_hydro_advance_extra_queue_c (void* packet, const int id) {
        return MILHOJA_ERROR_UNABLE_TO_RELEASE_STREAM;
        if (packet == nullptr) {
            std::cerr << "[release_DataPacket_gpu_tf_hydro_extra_queue_C] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        DataPacket_gpu_tf_hydro*   packet_h = static_cast<DataPacket_gpu_tf_hydro*>(packet);
    
        //try {
        //    packet_h->releaseExtraQueue(id);
        //} catch (const std::exception& exc) {
        //    std::cerr << exc.what() << std::endl;
        //    return MILHOJA_ERROR_UNABLE_TO_RELEASE_STREAM;
        //} catch (...) {
        //    std::cerr << "[release_DataPacket_gpu_tf_hydro_extra_queue_C] Unknown error caught" << std::endl;
        //    return MILHOJA_ERROR_UNABLE_TO_RELEASE_STREAM;
        //}
    
        return MILHOJA_SUCCESS;
    }

    //----- C TASK FUNCTION TO BE CALLED BY RUNTIME
    //------------------changed function name
    void dr_hydro_advance_packet_oacc_tf (const int threadIndex, void* dataItem_h) {
        DataPacket_gpu_tf_hydro* packet_h = static_cast<DataPacket_gpu_tf_hydro*>(dataItem_h);
        const int queue1_h = packet_h->asynchronousQueue();
        const int _nTiles_h = packet_h->_nTiles_h;
        // set queues to -1 for 2d problem
        const int queue2_h = -1;
        const int queue3_h = -1;
        
        //const int queue2_h = packet_h->extraAsynchronousQueue(2);
        //if (queue2_h < 0)
        //	throw std::overflow_error("[gpu_tf_hydro_Cpp2C] Potential overflow error when accessing async queue id.");
        //const int queue3_h = packet_h->extraAsynchronousQueue(3);
        //if (queue3_h < 0)
        //	throw std::overflow_error("[gpu_tf_hydro_Cpp2C] Potential overflow error when accessing async queue id.");
        

        void* _nTiles_d = static_cast<void*>( packet_h->_nTiles_d );
        void* _dt_d = static_cast<void*>( packet_h->_dt_d );
        void* _tile_deltas_d = static_cast<void*>( packet_h->_tile_deltas_d );
        void* _tile_lo_d = static_cast<void*>( packet_h->_tile_lo_d );
        void* _tile_hi_d = static_cast<void*>( packet_h->_tile_hi_d );
        void* _tile_lbound_d = static_cast<void*>( packet_h->_tile_lbound_d );
        void* _U_d = static_cast<void*>( packet_h->_U_d );
        void* _hydro_op1_auxc_d = static_cast<void*>( packet_h->_hydro_op1_auxc_d );
        void* _hydro_op1_flX_d = static_cast<void*>( packet_h->_hydro_op1_flX_d );
        void* _hydro_op1_flY_d = static_cast<void*>( packet_h->_hydro_op1_flY_d );
        void* _hydro_op1_flZ_d = static_cast<void*>( packet_h->_hydro_op1_flZ_d );
        

        // needed to change function name
        // Pass data packet info to C-to-Fortran Reinterpretation Layer
        dr_hydro_advance_packet_oacc_c2f (
        packet_h,
        queue1_h,
        queue2_h,
        queue3_h,
        _nTiles_h,
        _nTiles_d,
        _dt_d,
        _tile_deltas_d,
        _tile_lo_d,
        _tile_hi_d,
        _tile_lbound_d,
        _U_d,
        _hydro_op1_auxc_d,
        _hydro_op1_flX_d,
        _hydro_op1_flY_d,
        _hydro_op1_flZ_d
        
        );
    }
}
