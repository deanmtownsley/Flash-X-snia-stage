#include "Tile_cpu_tf_hydro.h"

#include <Milhoja_Runtime.h>
#include <Milhoja_RuntimeBackend.h>

void*  Tile_cpu_tf_hydro::scratch_hydro_op1_auxc_ = nullptr;
void*  Tile_cpu_tf_hydro::scratch_hydro_op1_flx_ = nullptr;
void*  Tile_cpu_tf_hydro::scratch_hydro_op1_fly_ = nullptr;
void*  Tile_cpu_tf_hydro::scratch_hydro_op1_flz_ = nullptr;

void Tile_cpu_tf_hydro::acquireScratch(void) {
    const unsigned int  nThreads = milhoja::Runtime::instance().nMaxThreadsPerTeam();

    if (scratch_hydro_op1_auxc_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::acquireScratch] scratch_hydro_op1_auxc_ scratch already allocated");
    } else if (scratch_hydro_op1_flx_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::acquireScratch] scratch_hydro_op1_flx_ scratch already allocated");
    } else if (scratch_hydro_op1_fly_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::acquireScratch] scratch_hydro_op1_fly_ scratch already allocated");
    } else if (scratch_hydro_op1_flz_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::acquireScratch] scratch_hydro_op1_flz_ scratch already allocated");
    }

    std::size_t nBytes = nThreads
             * Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_AUXC_SIZE_
             * sizeof(milhoja::Real);
    milhoja::RuntimeBackend::instance().requestCpuMemory(nBytes, &scratch_hydro_op1_auxc_);

    nBytes =   nThreads
             * Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_FLX_SIZE_
             * sizeof(milhoja::Real);
    milhoja::RuntimeBackend::instance().requestCpuMemory(nBytes, &scratch_hydro_op1_flx_);

    nBytes =   nThreads
             * Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_FLY_SIZE_
             * sizeof(milhoja::Real);
    milhoja::RuntimeBackend::instance().requestCpuMemory(nBytes, &scratch_hydro_op1_fly_);

    nBytes =   nThreads
             * Tile_cpu_tf_hydro::SCRATCH_HYDRO_OP1_FLZ_SIZE_
             * sizeof(milhoja::Real);
    milhoja::RuntimeBackend::instance().requestCpuMemory(nBytes, &scratch_hydro_op1_flz_);
}

void Tile_cpu_tf_hydro::releaseScratch(void) {
    if (!scratch_hydro_op1_auxc_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::releaseScratch] scratch_hydro_op1_auxc_ scratch not allocated");
    } else if (!scratch_hydro_op1_flx_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::releaseScratch] scratch_hydro_op1_flx_ scratch not allocated");
    } else if (!scratch_hydro_op1_fly_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::releaseScratch] scratch_hydro_op1_fly_ scratch not allocated");
    } else if (!scratch_hydro_op1_flz_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::releaseScratch] scratch_hydro_op1_flz_ scratch not allocated");
    }

    milhoja::RuntimeBackend::instance().releaseCpuMemory(&scratch_hydro_op1_auxc_);
    scratch_hydro_op1_auxc_ = nullptr;

    milhoja::RuntimeBackend::instance().releaseCpuMemory(&scratch_hydro_op1_flx_);
    scratch_hydro_op1_flx_ = nullptr;

    milhoja::RuntimeBackend::instance().releaseCpuMemory(&scratch_hydro_op1_fly_);
    scratch_hydro_op1_fly_ = nullptr;

    milhoja::RuntimeBackend::instance().releaseCpuMemory(&scratch_hydro_op1_flz_);
    scratch_hydro_op1_flz_ = nullptr;
}

Tile_cpu_tf_hydro::Tile_cpu_tf_hydro(const milhoja::Real external_hydro_op1_dt,
                                     const int external_hydro_op1_eosMode)
    : milhoja::TileWrapper{},
      external_hydro_op1_dt_{external_hydro_op1_dt},
      external_hydro_op1_eosMode_{external_hydro_op1_eosMode}
{
}

Tile_cpu_tf_hydro::~Tile_cpu_tf_hydro(void) {
}

std::unique_ptr<milhoja::TileWrapper> Tile_cpu_tf_hydro::clone(std::shared_ptr<milhoja::Tile>&& tileToWrap) const {
    Tile_cpu_tf_hydro* ptr = new Tile_cpu_tf_hydro{external_hydro_op1_dt_,
                                                   external_hydro_op1_eosMode_};

    if (ptr->tile_) {
        throw std::logic_error("[Tile_cpu_tf_hydro::clone] Internal tile_ member not null");
    }
    ptr->tile_ = std::move(tileToWrap);
    if (!(ptr->tile_) || tileToWrap) {
        throw std::logic_error("[Tile_cpu_tf_hydro::clone] Wrapper did not take ownership of tile");
    }

    return std::unique_ptr<milhoja::TileWrapper>{ptr};
}
