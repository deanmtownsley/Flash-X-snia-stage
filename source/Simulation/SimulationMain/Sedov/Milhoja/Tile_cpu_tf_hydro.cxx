#include "Tile_cpu_tf_hydro.h"

#include <Milhoja_Runtime.h>
#include <Milhoja_RuntimeBackend.h>
#ifdef DEBUG_RUNTIME
#include <Milhoja_Logger.h>
#endif

void*  Tile_cpu_tf_hydro::hydro_op1_auxc_ = nullptr;

void Tile_cpu_tf_hydro::acquireScratch(void) {
    const unsigned int  nThreads = milhoja::Runtime::instance().nMaxThreadsPerTeam();

    if (hydro_op1_auxc_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::acquireScratch] hydro_op1_auxc_ scratch already allocated");
    }

    const std::size_t nBytes = nThreads
                    * Tile_cpu_tf_hydro::HYDRO_OP1_AUXC_SIZE_
                    * sizeof(milhoja::Real);

    milhoja::RuntimeBackend::instance().requestCpuMemory(nBytes, &hydro_op1_auxc_);

#ifdef DEBUG_RUNTIME
    std::string   msg = "[Tile_cpu_tf_hydro::acquireScratch] Acquired"
                    + std::to_string(nThreads)
                    + " hydro_op1_auxc_ scratch blocks"
    milhoja::Logger::instance().log(msg);
#endif
}

void Tile_cpu_tf_hydro::releaseScratch(void) {
    if (!hydro_op1_auxc_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::releaseScratch] hydro_op1_auxc_ scratch not allocated");
    }

    milhoja::RuntimeBackend::instance().releaseCpuMemory(&hydro_op1_auxc_);
    hydro_op1_auxc_ = nullptr;

#ifdef DEBUG_RUNTIME
    std::string   msg = "[Tile_cpu_tf_hydro::releaseScratch] Released hydro_op1_auxc_ scratch"
    milhoja::Logger::instance().log(msg);
#endif
}

Tile_cpu_tf_hydro::Tile_cpu_tf_hydro(const milhoja::Real dt)
    : milhoja::TileWrapper{},
      dt_{dt}
{
#ifdef DEBUG_RUNTIME
    std::string   msg = "[Tile_cpu_tf_hydro] Creating wrapper object";
    milhoja::Logger::instance().log(msg);
#endif
}

Tile_cpu_tf_hydro::~Tile_cpu_tf_hydro(void) {
#ifdef DEBUG_RUNTIME
    std::string   msg = "[~Tile_cpu_tf_hydro] Destroying wrapper object";
    milhoja::Logger::instance().log(msg);
#endif
}

std::unique_ptr<milhoja::TileWrapper> Tile_cpu_tf_hydro::clone(std::shared_ptr<milhoja::Tile>&& tileToWrap) const {
    Tile_cpu_tf_hydro* ptr = new Tile_cpu_tf_hydro{dt_};

    if (ptr->tile_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::clone] Internal tile_ member not null");
    }
    ptr->tile_ = std::move(tileToWrap);
    if (!(ptr->tile_) || tileToWrap) {
        throw std::logic_error("[Tile_cpu_tf_hydro::clone] Wrapper did not take ownership of tile");
    }

    return std::unique_ptr<milhoja::TileWrapper>{ptr};
}
