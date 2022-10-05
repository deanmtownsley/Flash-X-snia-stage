#include "DataPacket_Hydro_gpu_3.h"

#include <cassert>
#include <cstring>
#include <stdexcept>

#include <Milhoja_IntVect.h>
#include <Milhoja_RealVect.h>
#include <Milhoja_Grid.h>
#include <Milhoja_RuntimeBackend.h>

/**
 * Construct a DataPacket containing no Tile objects and with no resources
 * assigned to it.
 *
 * @todo My gut feeling is that the contents of this file should be generated
 *       without knowledge of Sedov.h/Driver.h.  Can we get rid of those here?
 */
DataPacket_Hydro_gpu_3::DataPacket_Hydro_gpu_3(const milhoja::Real dt)
    : milhoja::DataPacket{},
#if MILHOJA_NDIM == 3
      stream2_{},
      stream3_{},
#endif
      nxb_{1},
      nyb_{1},
      nzb_{1},
      nCcVars_{milhoja::Grid::instance().getNCcVariables()},
      nFluxVars_{milhoja::Grid::instance().getNFluxVariables()},
      nGuard_{milhoja::Grid::instance().getNGuardcells()},
      nTiles_h_{0},
      nTiles_p_{nullptr},
      nTiles_d_{nullptr},
      dt_h_{dt},
      dt_p_{nullptr},
      dt_d_{nullptr},
      deltas_start_p_{nullptr},
      deltas_start_d_{nullptr},
      lo_start_p_{nullptr},
      lo_start_d_{nullptr},
      hi_start_p_{nullptr},
      hi_start_d_{nullptr},
      loGC_start_p_{nullptr},
      loGC_start_d_{nullptr},
      hiGC_start_p_{nullptr},
      hiGC_start_d_{nullptr},
      U_start_p_{nullptr},
      U_start_d_{nullptr},
      auxC_start_d_{nullptr},
      faceX_start_d_{nullptr},
      faceY_start_d_{nullptr},
      faceZ_start_d_{nullptr}
{
    // nyb_/nzb_ must be set to default value of 1 prior to this call so that
    // they have the correct value if NDIM < 3.
    milhoja::Grid::instance().getBlockSize(&nxb_, &nyb_, &nzb_);
}

/**
 * Destroy DataPacket.  Under normal circumstances, the DataPacket should have
 * been consumed and therefore own no resources.
 */
DataPacket_Hydro_gpu_3::~DataPacket_Hydro_gpu_3(void) {
#if MILHOJA_NDIM == 3
    if (stream2_.isValid() || stream3_.isValid()) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::~DataPacket_Hydro_gpu_3] "
                               "One or more extra streams not released");
    }
#endif
}

/**
 *
 */
std::unique_ptr<milhoja::DataPacket>   DataPacket_Hydro_gpu_3::clone(void) const {
    return std::unique_ptr<milhoja::DataPacket>{new DataPacket_Hydro_gpu_3{dt_h_}};
}

/**
 *
 * @todo Error check pointers
 */
void  DataPacket_Hydro_gpu_3::tileSize_host(int* nxbGC,
                                            int* nybGC,
                                            int* nzbGC,
                                            int* nCcVars,
                                            int* nFluxVars) const {
    *nxbGC = static_cast<int>(nxb_ + 2 * nGuard_ * MILHOJA_K1D);
    *nybGC = static_cast<int>(nyb_ + 2 * nGuard_ * MILHOJA_K2D);
    *nzbGC = static_cast<int>(nzb_ + 2 * nGuard_ * MILHOJA_K3D);
    // We are not including GAME in U, which is the last variable in each block
    // Exclude all the analytic variables as well, which should be after GAME.
    *nCcVars = nCcVars_ - 8;
    *nFluxVars = nFluxVars_;
}

#if MILHOJA_NDIM == 3 && defined(MILHOJA_OPENACC_OFFLOADING)
/**
 * Refer to the documentation of this member function for DataPacket.
 */
void  DataPacket_Hydro_gpu_3::releaseExtraQueue(const unsigned int id) {
    if        (id == 2) {
        if (!stream2_.isValid()) {
            throw std::logic_error("[DataPacket_Hydro_gpu_3::releaseExtraQueue] "
                                   "Second queue invalid or already released");
        } else {
            milhoja::RuntimeBackend::instance().releaseStream(stream2_);
        }
    } else if (id == 3) {
        if (!stream3_.isValid()) {
            throw std::logic_error("[DataPacket_Hydro_gpu_3::releaseExtraQueue] "
                                   "Third queue invalid or already released");
        } else {
            milhoja::RuntimeBackend::instance().releaseStream(stream3_);
        }
    } else {
        throw std::invalid_argument("[DataPacket_Hydro_gpu_3::releaseExtraQueue] "
                                    "Invalid id");
    }
}
#endif

#if MILHOJA_NDIM == 3 && defined(MILHOJA_OPENACC_OFFLOADING)
/**
 * Refer to the documentation of this member function for DataPacket.
 */
int  DataPacket_Hydro_gpu_3::extraAsynchronousQueue(const unsigned int id) {
    if        (id == 2) {
        if (!stream2_.isValid()) {
            throw std::logic_error("[DataPacket_Hydro_gpu_3::extraAsynchronousQueue] "
                                   "Second queue invalid");
        } else {
            return stream2_.accAsyncQueue;
        }
    } else if (id == 3) {
        if (!stream3_.isValid()) {
            throw std::logic_error("[DataPacket_Hydro_gpu_3::extraAsynchronousQueue] "
                                   "Third queue invalid");
        } else {
            return stream3_.accAsyncQueue;
        }
    } else {
        throw std::invalid_argument("[DataPacket_Hydro_gpu_3::extraAsynchronousQueue] Invalid id");
    }
}
#endif

/**
 *
 * @todo How do variables such as dt get into this code?
 */
void  DataPacket_Hydro_gpu_3::pack(void) {
    using namespace milhoja;

    std::string   errMsg = isNull();
    if (errMsg != "") {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::pack] " + errMsg);
    } else if (tiles_.size() == 0) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::pack] No tiles added");
    }

    //----- OBTAIN NON-TILE-SPECIFIC HOST-SIDE DATA
    // TODO: Check for correctness of cast here of elsewhere?
    // This cast is necessary since Fortran code consumes the packet.
    nTiles_h_ = static_cast<int>(tiles_.size());
    int   nxbGC_h     = -1;
    int   nybGC_h     = -1;
    int   nzbGC_h     = -1;
    int   nCcVars_h   = -1;
    int   nFluxVars_h = -1;
    tileSize_host(&nxbGC_h, &nybGC_h, &nzbGC_h, &nCcVars_h, &nFluxVars_h);

    //----- COMPUTE SIZES/OFFSETS & ACQUIRE MEMORY
    // TODO: Johann would like to minimize the use of pinned memory.  Presently,
    // we structure a packet in pinned memory, send the data to the device, and
    // copy data back in accord with that same structure.  There is no need to
    // do this.  Rather, we can copy back data to pinned host with its own
    // structure and overwrite the previous structure.
    constexpr std::size_t    sz_nTiles = sizeof(int);
    constexpr std::size_t    sz_dt     = sizeof(Real);

    constexpr std::size_t    sz_deltas = MILHOJA_MDIM * sizeof(Real);
    constexpr std::size_t    sz_lo     = MILHOJA_MDIM * sizeof(int);
    constexpr std::size_t    sz_hi     = MILHOJA_MDIM * sizeof(int);
    constexpr std::size_t    sz_loGC   = MILHOJA_MDIM * sizeof(int);
    constexpr std::size_t    sz_hiGC   = MILHOJA_MDIM * sizeof(int);
    std::size_t              sz_U      =  nxbGC_h * nybGC_h * nzbGC_h * nCcVars_h * sizeof(Real);
    std::size_t              sz_auxC   =  nxbGC_h * nybGC_h * nzbGC_h *             sizeof(Real);
    // TODO: Can't this be sized smaller for simpleUnsplit?
    std::size_t              sz_FX     = (nxbGC_h + 1) *  nybGC_h      *  nzbGC_h      * nFluxVars_h * sizeof(Real);
#if MILHOJA_NDIM == 1
    std::size_t              sz_FY     =  1            *  1            *  1            * 1           * sizeof(Real);
#else
    std::size_t              sz_FY     =  nxbGC_h      * (nybGC_h + 1) *  nzbGC_h      * nFluxVars_h * sizeof(Real);
#endif
#if MILHOJA_NDIM <= 2
    // The update solution function is passed FZ regardless of the fact that it
    // won't use it.  As part of this, we get FZ for the nth tile.  Therefore,
    // we have a single element per tile so that indexing into FZ by tile still
    // yields a valid array.  We could make this array a single element, but the
    // code generator would need to write that code so that it only every
    // indexed into FZ for the first tile.
    // TODO: Is this still necessary?
    std::size_t              sz_FZ     =  1            *  1            *  1            * 1           * sizeof(Real);
#else
    std::size_t              sz_FZ     =  nxbGC_h      *  nybGC_h      * (nzbGC_h + 1) * nFluxVars_h * sizeof(Real);
#endif

    std::size_t    nScratchBytes = nTiles_h_ * (sz_auxC + sz_FX + sz_FY + sz_FZ);
    std::size_t    nScratchBytes_padded = pad(nScratchBytes);
    if ((nScratchBytes_padded % ALIGN_SIZE) != 0) {
        throw std::logic_error("[Packet] nScratchBytes padding failed");
    }

    std::size_t    nCopyInBytes =   sz_nTiles
                                  + sz_dt
                                  + nTiles_h_
                                    * (  sz_deltas
                                       + sz_lo   + sz_hi
                                       + sz_loGC + sz_hiGC);
    std::size_t    nCopyInBytes_padded = pad(nCopyInBytes);
    if ((nCopyInBytes_padded % ALIGN_SIZE) != 0) {
        throw std::logic_error("[Packet] nCopyInBytes padding failed");
    }

    std::size_t    nCopyInOutBytes = nTiles_h_ * sz_U;
    std::size_t    nCopyInOutBytes_padded = pad(nCopyInOutBytes);
    if ((nCopyInOutBytes_padded % ALIGN_SIZE) != 0) {
        throw std::logic_error("[Packet] nCopyInOutBytes padding failed");
    }

    std::size_t    nCopyOutBytes   = 0;
    std::size_t    nCopyOutBytes_padded = pad(nCopyOutBytes);
    if ((nCopyOutBytes_padded % ALIGN_SIZE) != 0) {
        throw std::logic_error("[Packet] nCopyOutBytes padding failed");
    }

    nCopyToGpuBytes_    = nCopyInBytes_padded    + nCopyInOutBytes;
    nReturnToHostBytes_ = nCopyInOutBytes_padded + nCopyOutBytes;
    std::size_t  nBytesPerPacket =   nScratchBytes_padded
                                   + nCopyInBytes_padded
                                   + nCopyInOutBytes_padded
                                   + nCopyOutBytes_padded;

    // ACQUIRE PINNED AND GPU MEMORY & SPECIFY STRUCTURE
    // Scratch only needed on GPU side
    // At present, this call makes certain that each memory buffer acquires is
    // appropriately byte aligned.
    RuntimeBackend::instance().requestGpuMemory(nBytesPerPacket - nScratchBytes_padded,
                                                &packet_p_, nBytesPerPacket, &packet_d_);

    //----- BYTE-ALIGN SCRATCH SECTION
    // Order from largest to smallest in data type size
    static_assert(sizeof(char) == 1, "Invalid char size");
    char*   ptr_d = static_cast<char*>(packet_d_);

    auxC_start_d_  = static_cast<void*>(ptr_d);
    ptr_d += nTiles_h_ * sz_auxC;

    faceX_start_d_ = static_cast<void*>(ptr_d);
    ptr_d += nTiles_h_ * sz_FX;
    faceY_start_d_ = static_cast<void*>(ptr_d);
    ptr_d += nTiles_h_ * sz_FY;
    faceZ_start_d_ = static_cast<void*>(ptr_d);
    ptr_d += nTiles_h_ * sz_FZ;
 
    //----- BYTE-ALIGN COPY-IN SECTION
    // Order from largest to smallest in data type size
    copyInStart_p_ =   static_cast<char*>(packet_p_);
    copyInStart_d_ =   static_cast<char*>(packet_d_)
                     + nScratchBytes_padded;

    char*   ptr_p = copyInStart_p_;
    ptr_d         = copyInStart_d_;

    //-- REALS (8-byte)
    // non-tile-specific reals
    dt_p_ = static_cast<void*>(ptr_p);
    dt_d_ = static_cast<void*>(ptr_d);
    ptr_p += sz_dt;
    ptr_d += sz_dt;

    // tile-specific reals
    deltas_start_p_ = static_cast<void*>(ptr_p);
    deltas_start_d_ = static_cast<void*>(ptr_d);
    ptr_p += nTiles_h_ * sz_deltas;
    ptr_d += nTiles_h_ * sz_deltas;

    //-- INTEGERS (4-byte)
    // non-tile-specific integers
    nTiles_p_ = static_cast<void*>(ptr_p);
    nTiles_d_ = static_cast<void*>(ptr_d);
    ptr_p += sz_nTiles;
    ptr_d += sz_nTiles;

    // tile-specific integers
    lo_start_p_ = static_cast<void*>(ptr_p);
    lo_start_d_ = static_cast<void*>(ptr_d);
    ptr_p += nTiles_h_ * sz_lo;
    ptr_d += nTiles_h_ * sz_lo;

    hi_start_p_ = static_cast<void*>(ptr_p);
    hi_start_d_ = static_cast<void*>(ptr_d);
    ptr_p += nTiles_h_ * sz_hi;
    ptr_d += nTiles_h_ * sz_hi;

    loGC_start_p_ = static_cast<void*>(ptr_p);
    loGC_start_d_ = static_cast<void*>(ptr_d);
    ptr_p += nTiles_h_ * sz_loGC;
    ptr_d += nTiles_h_ * sz_loGC;

    hiGC_start_p_ = static_cast<void*>(ptr_p);
    hiGC_start_d_ = static_cast<void*>(ptr_d);
    ptr_p += nTiles_h_ * sz_hiGC;
    ptr_d += nTiles_h_ * sz_hiGC;

    //----- BYTE-ALIGN COPY-IN/-OUT SECTION
    // Order from largest to smallest in data type size
    // 
    // Pad copy-in section if necessary to get correct byte alignment
    location_ = PacketDataLocation::CC1;

    copyInOutStart_p_ = copyInStart_p_ + nCopyInBytes_padded;
    copyInOutStart_d_ = copyInStart_d_ + nCopyInBytes_padded;

    ptr_p = static_cast<char*>(copyInOutStart_p_);
    ptr_d = static_cast<char*>(copyInOutStart_d_);

    U_start_p_ = static_cast<void*>(ptr_p);
    U_start_d_ = static_cast<void*>(ptr_d);
    ptr_p += nTiles_h_ * sz_U;
    ptr_d += nTiles_h_ * sz_U;

    // No copy-out data

    // Store for later unpacking the location in pinned memory of the different
    // blocks.
    if (pinnedPtrs_) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::pack] Pinned pointers already exist");
    }
    pinnedPtrs_ = new BlockPointersPinned[nTiles_h_];

    //----- SCRATCH SECTION
    // Nothing to include nor record

    //----- COPY IN SECTION
    // Non-tile-specific data
    std::memcpy(nTiles_p_, static_cast<void*>(&nTiles_h_), sz_nTiles);
    std::memcpy(dt_p_,     static_cast<void*>(&dt_h_),     sz_dt);

    // Tile-specific metadata
    char*   char_ptr;
    for (std::size_t n=0; n<nTiles_h_; ++n) {
        Tile*   tileDesc_h = tiles_[n].get();
        if (tileDesc_h == nullptr) {
            throw std::runtime_error("[DataPacket_Hydro_gpu_3::pack] Bad tileDesc");
        }

        const RealVect      deltas = tileDesc_h->deltas();
        const IntVect       lo     = tileDesc_h->lo();
        const IntVect       hi     = tileDesc_h->hi();
        const IntVect       loGC   = tileDesc_h->loGC();
        const IntVect       hiGC   = tileDesc_h->hiGC();
        Real*               data_h = tileDesc_h->dataPtr();
        if (data_h == nullptr) {
            throw std::logic_error("[DataPacket_Hydro_gpu_3::pack] "
                                   "Invalid pointer to data in host memory");
        }

        Real    deltas_h[MILHOJA_MDIM] = {deltas.I(), deltas.J(), deltas.K()};
        char_ptr = static_cast<char*>(deltas_start_p_) + n * sz_deltas;
        std::memcpy(static_cast<void*>(char_ptr),
                    static_cast<void*>(deltas_h),
                    sz_deltas);

        // Global index space is 0-based in runtime; 1-based in Fortran code.
        // Translate here so that it is immediately ready for use with Fortran.
        int     lo_h[MILHOJA_MDIM] = {lo.I()+1, lo.J()+1, lo.K()+1};
        char_ptr = static_cast<char*>(lo_start_p_) + n * sz_lo;
        std::memcpy(static_cast<void*>(char_ptr),
                    static_cast<void*>(lo_h),
                    sz_lo);

        // Global index space is 0-based in runtime; 1-based in Fortran code.
        // Translate here so that it is immediately ready for use with Fortran.
        int     hi_h[MILHOJA_MDIM] = {hi.I()+1, hi.J()+1, hi.K()+1};
        char_ptr = static_cast<char*>(hi_start_p_) + n * sz_hi;
        std::memcpy(static_cast<void*>(char_ptr),
                    static_cast<void*>(hi_h),
                    sz_hi);

        // Global index space is 0-based in runtime; 1-based in Fortran code.
        // Translate here so that it is immediately ready for use with Fortran.
        int     loGC_h[MILHOJA_MDIM] = {loGC.I()+1, loGC.J()+1, loGC.K()+1};
        char_ptr = static_cast<char*>(loGC_start_p_) + n * sz_loGC;
        std::memcpy(static_cast<void*>(char_ptr),
                    static_cast<void*>(loGC_h),
                    sz_loGC);

        // Global index space is 0-based in runtime; 1-based in Fortran code.
        // Translate here so that it is immediately ready for use with Fortran.
        int     hiGC_h[MILHOJA_MDIM] = {hiGC.I()+1, hiGC.J()+1, hiGC.K()+1};
        char_ptr = static_cast<char*>(hiGC_start_p_) + n * sz_hiGC;
        std::memcpy(static_cast<void*>(char_ptr),
                    static_cast<void*>(hiGC_h),
                    sz_hiGC);

        char_ptr = static_cast<char*>(U_start_p_) + n * sz_U;
        std::memcpy(static_cast<void*>(char_ptr),
                    static_cast<void*>(data_h),
                    sz_U);
        pinnedPtrs_[n].CC1_data = static_cast<Real*>((void*)char_ptr);
        // Data will always be copied back from CC1
        pinnedPtrs_[n].CC2_data = nullptr;
    }

    // Request memory first and pack data in pinned memory *before* acquiring
    // streams.  This is a useful optimization if stream resources are too
    // little and the pack() routine gets blocked waiting for stream resources.
    // With this ordering, we pack the memory, which can be slow and have the
    // packet ready once a stream is available.  Indeed, in some cases a stream
    // might not be ready before we pack, but will be ready after we pack.
    stream_ = RuntimeBackend::instance().requestStream(true);
    if (!stream_.isValid()) {
        throw std::runtime_error("[DataPacket_Hydro_gpu_3::pack] Unable to acquire stream");
    }
#if MILHOJA_NDIM == 3
    stream2_ = RuntimeBackend::instance().requestStream(true);
    stream3_ = RuntimeBackend::instance().requestStream(true);
    if (!stream2_.isValid() || !stream3_.isValid()) {
        throw std::runtime_error("[DataPacket_Hydro_gpu_3::pack] Unable to acquire extra streams");
    }
#endif
}

/**
 * The runtime calls this member function automatically once a DataPacket
 * has arrived in the host memory again.  It is responsible for unpacking
 * the contents and in particular for copying cell-centered data back to the
 * host-side Grid data structures that hold solution data.  The data is copied
 * back in accord with the variable masks set in the DataPacket to avoid
 * inadvertently overwriting variables that were updated in parallel by other
 * actions.
 *
 * All memory and stream resources are released.
 *
 * While the packet is consumed once the function finishes, the list of Tiles
 * that were included in the packet is preserved.  This is necessary so that
 * runtime elements such as MoverUnpacker can enqueue the Tiles with its data
 * subscriber.
 *
 * @todo Should unpacking be made more generic so that the CC blocks need not 
 *       start always with the first data variable.  What if the packet just
 *       needs to include variables 3-5 (out of 10 for example)?
 */
void  DataPacket_Hydro_gpu_3::unpack(void) {
    using namespace milhoja;

    if (tiles_.size() <= 0) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] "
                               "Empty data packet");
    } else if (!stream_.isValid()) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] "
                               "Stream not acquired");
    } else if (pinnedPtrs_ == nullptr) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] "
                               "No pinned pointers set");
    } else if (   (startVariable_ < 0) || (startVariable_ >= nCcVars_)
               || (endVariable_   < 0) || (endVariable_   >= nCcVars_)) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] "
                               "Invalid variable mask");
    } else if (startVariable_ > endVariable_) {
        throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] "
                               "start variable ID > end variable ID");
    }

    // Release stream as soon as possible
    RuntimeBackend::instance().releaseStream(stream_);
    assert(!stream_.isValid());

    unsigned int   N_ELEMENTS_PER_CC_PER_VARIABLE =   (nxb_ + 2 * nGuard_ * MILHOJA_K1D)
                                                    * (nyb_ + 2 * nGuard_ * MILHOJA_K2D)
                                                    * (nzb_ + 2 * nGuard_ * MILHOJA_K3D);

    for (auto n=0; n<tiles_.size(); ++n) {
        Tile*   tileDesc_h = tiles_[n].get();

        Real*         data_h = tileDesc_h->dataPtr();
        const Real*   data_p = nullptr;
        switch (location_) {
            case PacketDataLocation::CC1:
                data_p = pinnedPtrs_[n].CC1_data;
                break;
            case PacketDataLocation::CC2:
                data_p = pinnedPtrs_[n].CC2_data;
                break;
            default:
                throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] Data not in CC1 or CC2");
        }

        if (data_h == nullptr) {
            throw std::logic_error("[DataPacket_Hydro_gpu_3::unpack] "
                                   "Invalid pointer to data in host memory");
        } else if (data_p == nullptr) {
            throw std::runtime_error("[DataPacket_Hydro_gpu_3::unpack] "
                                     "Invalid pointer to data in pinned memory");
        }

        // The code here imposes requirements on the variable indices.
        std::size_t  offset =   N_ELEMENTS_PER_CC_PER_VARIABLE
                              * static_cast<std::size_t>(startVariable_);
        Real*        start_h = data_h + offset;
        const Real*  start_p = data_p + offset;
        std::size_t  nBytes =  (endVariable_ - startVariable_ + 1)
                              * N_ELEMENTS_PER_CC_PER_VARIABLE
                              * sizeof(Real);
        std::memcpy((void*)start_h, (void*)start_p, nBytes);
    }
}

