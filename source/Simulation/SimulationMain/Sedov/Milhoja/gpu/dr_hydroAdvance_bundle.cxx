#include <Milhoja.h>
#include <Milhoja_real.h>
#include <Milhoja_interface_error_codes.h>

#include "DataPacket_Hydro_gpu_3.h"

// TODO: Code generated for use by the runtime should not depend on the
// application at all.  Therefore, the need to include this header is
// unacceptable.  How to get variable mask information otherwise?
#include "Simulation.h"

#ifndef MILHOJA_OPENACC_OFFLOADING
#error "This file should only be compiled if using OpenACC offloading"
#endif

extern "C" {
    //----- C DECLARATION OF FORTRAN ROUTINE WITH C-COMPATIBLE INTERFACE
    void   dr_hydro_advance_packet_oacc_c2f(void* packet_h,
                                            const int dataQ_h,
                                            const int queue2_h,
                                            const int queue3_h,
                                            const int nTiles_h,
                                            const int nxbGC_h, const int nybGC_h, const int nzbGC_h,
                                            const int nCcVar_h, const int nFluxVar_h,
                                            const void* nTiles_d, const void* dt_d,
                                            const void* deltas_start_d,
                                            const void* lo_start_d,   const void* hi_start_d,
                                            const void* loGC_start_d,
                                            const void* U_start_d,
                                            const void* auxC_start_d,
                                            const void* faceX_start_d,
                                            const void* faceY_start_d,
                                            const void* faceZ_start_d);

    int instantiate_hydro_advance_packet_c(const milhoja::Real dt, void** packet) {
        if        ( packet == nullptr) {
            std::cerr << "[instantiate_hydro_advance_packet_c] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        } else if (*packet != nullptr) {
            std::cerr << "[instantiate_hydro_advance_packet_c] *packet not NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_NOT_NULL;
        }

        try {
            *packet = static_cast<void*>(new DataPacket_Hydro_gpu_3(dt));
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_PACKET;
        } catch (...) {
            std::cerr << "[instantiate_hydro_advance_packet_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_CREATE_PACKET;
        }

        return MILHOJA_SUCCESS;
    }

    // TODO: Should this take a pointer to a pointer so that packet can be
    // nulled by this function?
    int delete_hydro_advance_packet_c(void* packet) {
        if (packet == nullptr) {
            std::cerr << "[delete_hydro_advance_packet_c] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        delete static_cast<DataPacket_Hydro_gpu_3*>(packet);

        return MILHOJA_SUCCESS;
    }

    int release_hydro_advance_extra_queue_c(void* packet, const int id) {
        if (packet == nullptr) {
            std::cerr << "[release_hydro_advance_extra_queue_c] packet is NULL" << std::endl;
            return MILHOJA_ERROR_POINTER_IS_NULL;
        }
        DataPacket_Hydro_gpu_3*   packet_h = static_cast<DataPacket_Hydro_gpu_3*>(packet);

        try {
            packet_h->releaseExtraQueue(id);
        } catch (const std::exception& exc) {
            std::cerr << exc.what() << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_RELEASE_STREAM;
        } catch (...) {
            std::cerr << "[release_hydro_advance_extra_queue_c] Unknown error caught" << std::endl;
            return MILHOJA_ERROR_UNABLE_TO_RELEASE_STREAM;
        }

        return MILHOJA_SUCCESS;
    }

    //----- C TASK FUNCTION TO BE CALLED BY RUNTIME
    void dr_hydro_advance_packet_oacc_tf(const int tId, void* dataItem_h) {
        DataPacket_Hydro_gpu_3*             packet_h = static_cast<DataPacket_Hydro_gpu_3*>(dataItem_h);
        const milhoja::PacketDataLocation   location = packet_h->getDataLocation();
        const int                           dataQ_h  = packet_h->asynchronousQueue();
#if MILHOJA_NDIM == 3
        const int                           queue2_h = packet_h->extraAsynchronousQueue(2);
        const int                           queue3_h = packet_h->extraAsynchronousQueue(3);
#else
        const int                           queue2_h = -1;
        const int                           queue3_h = -1;
#endif
        const int                           nTiles_h = packet_h->nTiles_host();
        int   nxbGC_h  = -1;
        int   nybGC_h  = -1;
        int   nzbGC_h  = -1;
        int   nCcVar_h = -1;
        int   nFluxVar_h = -1;
        packet_h->tileSize_host(&nxbGC_h, &nybGC_h, &nzbGC_h, &nCcVar_h, &nFluxVar_h);

        int*              nTiles_d       = packet_h->nTiles_devptr();
        milhoja::Real*    dt_d           = packet_h->dt_devptr();
        milhoja::Real*    deltas_start_d = packet_h->deltas_devptr();
        int*              lo_start_d     = packet_h->lo_devptr();
        int*              hi_start_d     = packet_h->hi_devptr();
        int*              loGC_start_d   = packet_h->loGC_devptr();
        milhoja::Real*    U_start_d      = packet_h->U_devptr();
        milhoja::Real*    auxC_start_d   = packet_h->scratchAuxC_devptr();
        milhoja::Real*    faceX_start_d  = packet_h->scratchFaceX_devptr();
        milhoja::Real*    faceY_start_d  = packet_h->scratchFaceY_devptr();
        milhoja::Real*    faceZ_start_d  = packet_h->scratchFaceZ_devptr();

        // This task function neither reads from nor writes to GAME.  While it does
        // read from GAMC, this variable is not written to as part of the task
        // function's work.  Therefore, GAME need not be included in the packet and
        // GAMC need not be copied back to Grid data structures as part of
        // host-side unpacking.
        // 
        // For some versions of this task function, the following masking of
        // variables is not an optimization.  Without this masking, whatever
        // data was originally in GAMC/GAME in a packet might be used to
        // overwrite true values for these two variables during host-side
        // unpacking.  
        //
        // Note that to avoid such overwriting, GAMC must be adjacent in memory
        // to all other variables in the packet and GAME outside of this grouping.
        // For this test, these two variables were declared in Simulation.h as the
        // last two UNK variables to accomplish this goal.
        //
        // TODO: How to do the masking?  Does the setup tool/offline toolchain have
        // to determine how to assign indices to the variables so that this can
        // happen for all task functions that must filter?  Selecting the order of
        // variables in memory sounds like part of the larger optimization problem
        // as it affects all data packets.
        packet_h->setVariableMask(UNK_VARS_BEGIN-1, EINT_VAR-1);

        if (location != milhoja::PacketDataLocation::CC1) {
            throw std::runtime_error("[hy_advance_solution_packet_oacc_tf] "
                                     "Input data must be in CC1");
        }

        // Pass data packet info to C-to-Fortran Reinterpretation Layer
        dr_hydro_advance_packet_oacc_c2f(packet_h,
                                         dataQ_h, queue2_h, queue3_h,
                                         nTiles_h,
                                         nxbGC_h, nybGC_h, nzbGC_h,
                                         nCcVar_h, nFluxVar_h,
                                         nTiles_d, dt_d,
                                         deltas_start_d,
                                         lo_start_d,   hi_start_d,
                                         loGC_start_d,
                                         U_start_d,
                                         auxC_start_d,
                                         faceX_start_d, faceY_start_d, faceZ_start_d);
    }
}

