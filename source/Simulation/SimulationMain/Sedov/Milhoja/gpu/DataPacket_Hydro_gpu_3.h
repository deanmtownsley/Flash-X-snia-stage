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

#ifndef DATA_PACKET_HYDRO_GPU_3_H__
#define DATA_PACKET_HYDRO_GPU_3_H__

#include <Milhoja.h>
#include <Milhoja_real.h>
#include <Milhoja_DataPacket.h> 

class DataPacket_Hydro_gpu_3 : public milhoja::DataPacket {
public:
    std::unique_ptr<milhoja::DataPacket>  clone(void) const override;

    DataPacket_Hydro_gpu_3(const milhoja::Real dt);
    ~DataPacket_Hydro_gpu_3(void);

    DataPacket_Hydro_gpu_3(DataPacket_Hydro_gpu_3&)                  = delete;
    DataPacket_Hydro_gpu_3(const DataPacket_Hydro_gpu_3&)            = delete;
    DataPacket_Hydro_gpu_3(DataPacket_Hydro_gpu_3&& packet)          = delete;
    DataPacket_Hydro_gpu_3& operator=(DataPacket_Hydro_gpu_3&)       = delete;
    DataPacket_Hydro_gpu_3& operator=(const DataPacket_Hydro_gpu_3&) = delete;
    DataPacket_Hydro_gpu_3& operator=(DataPacket_Hydro_gpu_3&& rhs)  = delete;

    void    pack(void) override;
    void    unpack(void) override;

#if MILHOJA_NDIM == 3 && defined(MILHOJA_OPENACC_OFFLOADING)
    int     extraAsynchronousQueue(const unsigned int id) override;
    void    releaseExtraQueue(const unsigned int id) override;
#endif

    // Since this is written for Fortran, might as well use int in the interface
    // rather than unsigned int.
    int              nTiles_host(void) const
        { return nTiles_h_; }
    milhoja::Real    dt_host(void) const
        { return dt_h_; }
    void             tileSize_host(int* nxbGC_h,
                                   int* nybGC_h,
                                   int* nzbGC_h,
                                   int* nCCVars_h,
                                   int* nFluxVars_h) const;

    int*             nTiles_devptr(void) const
        { return static_cast<int*>(nTiles_d_); }
    milhoja::Real*   dt_devptr(void) const
        { return static_cast<milhoja::Real*>(dt_d_); }

    milhoja::Real*   deltas_devptr(void) const
        { return static_cast<milhoja::Real*>(deltas_start_d_); }
    int*             lo_devptr(void) const
        { return static_cast<int*>(lo_start_d_); }
    int*             hi_devptr(void) const
        { return static_cast<int*>(hi_start_d_); }
    int*             loGC_devptr(void) const
        { return static_cast<int*>(loGC_start_d_); }
    int*             hiGC_devptr(void) const
        { return static_cast<int*>(hiGC_start_d_); }
    milhoja::Real*   U_devptr(void) const
        { return static_cast<milhoja::Real*>(U_start_d_); }
    milhoja::Real*   scratchAuxC_devptr(void) const
        { return static_cast<milhoja::Real*>(auxC_start_d_); }
    milhoja::Real*   scratchFaceX_devptr(void) const
        { return static_cast<milhoja::Real*>(faceX_start_d_); }
    milhoja::Real*   scratchFaceY_devptr(void) const
        { return static_cast<milhoja::Real*>(faceY_start_d_); }
    milhoja::Real*   scratchFaceZ_devptr(void) const
        { return static_cast<milhoja::Real*>(faceZ_start_d_); }

private:
    // Specify byte alignment of each memory request
    static constexpr std::size_t    ALIGN_SIZE = 16;
    static constexpr std::size_t    pad(const std::size_t sz) 
        { return ((sz + ALIGN_SIZE - 1) / ALIGN_SIZE) * ALIGN_SIZE; }

#if MILHOJA_NDIM == 3 && defined(MILHOJA_OPENACC_OFFLOADING)
    milhoja::Stream  stream2_;
    milhoja::Stream  stream3_;
#endif

    unsigned int           nxb_;
    unsigned int           nyb_;
    unsigned int           nzb_;
    const unsigned int     nCcVars_;
    const unsigned int     nFluxVars_;
    const unsigned int     nGuard_;

    int              nTiles_h_;
    void*            nTiles_p_;
    void*            nTiles_d_;

    milhoja::Real    dt_h_;
    void*            dt_p_;
    void*            dt_d_;

    void*            deltas_start_p_;
    void*            deltas_start_d_;

    void*            lo_start_p_;
    void*            lo_start_d_;
    void*            hi_start_p_;
    void*            hi_start_d_;

    void*            loGC_start_p_;
    void*            loGC_start_d_;
    void*            hiGC_start_p_;
    void*            hiGC_start_d_;

    void*            U_start_p_;
    void*            U_start_d_;

    void*            auxC_start_d_;
    void*            faceX_start_d_;
    void*            faceY_start_d_;
    void*            faceZ_start_d_;
};

#endif

