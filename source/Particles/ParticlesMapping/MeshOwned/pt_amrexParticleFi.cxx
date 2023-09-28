// Adapted from https://github.com/AMReX-Codes/amrex/blob/development/Src/F_Interfaces/Particle/AMReX_particlecontainer_fi.cpp
#include <AMReX_AmrParticles.H>
#include <AMReX_AmrCore.H>
#include <AMReX_ParallelDescriptor.H>
#include "constants.h"
#include "Simulation.h"

using namespace amrex;

#ifndef AMREX_FI_NSTRUCTREAL

#ifdef TYPE_PART_PROP
#define AMREX_FI_NSTRUCTREAL (BL_SPACEDIM+NPART_PROPS-(2*MDIM+5))
#else
#define AMREX_FI_NSTRUCTREAL (BL_SPACEDIM+NPART_PROPS-(2*MDIM+4))
#endif
#endif

#define AMREX_FI_NSTRUCTINT 0

namespace {
    using FParticleContainer = AmrParticleContainer<AMREX_FI_NSTRUCTREAL,
                                                    AMREX_FI_NSTRUCTINT>;
}

extern "C" {

    void pt_amrex_fi_new_particlecontainer (FParticleContainer*& particlecontainer,
                                         AmrCore* amrcore)
    {
	particlecontainer = new FParticleContainer(amrcore);
    }

    void pt_amrex_fi_delete_particlecontainer (FParticleContainer* particlecontainer)
    {
	delete particlecontainer;
    }

    void pt_amrex_fi_get_next_particle_id (Long& id)
    {
        id = FParticleContainer::ParticleType::NextID();
    }

    void pt_amrex_fi_get_cpu (int& cpu)
    {
        cpu = ParallelDescriptor::MyProc();
    }

    void pt_amrex_fi_get_particle_id(Long& id, const FParticleContainer::ParticleType* p)
    {
        id = p->id();
    }

    void pt_amrex_fi_set_particle_id(const Long& id, FParticleContainer::ParticleType* p)
    {
        p->id() = id;
    }

    void pt_amrex_fi_get_particle_cpu(int& cpu, const FParticleContainer::ParticleType* p)
    {
        cpu = p->cpu();
    }

    void pt_amrex_fi_set_particle_cpu(const int& cpu, FParticleContainer::ParticleType* p)
    {
        p->cpu() = cpu;
    }

    void pt_amrex_fi_write_particles(FParticleContainer* particlecontainer,
                                  const char* dirname, const char* pname, int is_checkpoint)
    {
        particlecontainer->Checkpoint(dirname, pname, is_checkpoint);
    }

    void pt_amrex_fi_particle_redistribute (FParticleContainer* particlecontainer,
                                         int lev_min, int lev_max, int ng)
    {
	particlecontainer->Redistribute(lev_min, lev_max, ng);
    }

    void pt_amrex_fi_get_particles_mfi(FParticleContainer* particlecontainer,
                                    int lev, MFIter* mfi, ParticleReal*& dp, Long& np)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
            if (np > 0) {
                auto& aos = particle_tile.GetArrayOfStructs();
                dp = reinterpret_cast<ParticleReal*>(aos.data());
            } else {
                dp = nullptr;
            }
        } else {
            np = 0;
            dp = nullptr;
        }
    }

    void pt_amrex_fi_add_particle_mfi(FParticleContainer* particlecontainer,
                                   int lev, MFIter* mfi, FParticleContainer::ParticleType* p)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto& particle_tile  = particle_level[std::make_pair(grid, tile)];
        particle_tile.push_back(*p);
    }

    void pt_amrex_fi_num_particles_mfi(FParticleContainer* particlecontainer,
                                    int lev, MFIter* mfi, Long& np)
    {
        const int grid = mfi->index();
        const int tile = mfi->LocalTileIndex();
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
        } else {
            np = 0;
        }
    }

    void pt_amrex_fi_get_particles_i(FParticleContainer* particlecontainer,
                                  int lev, int grid, int tile, ParticleReal*& dp, Long& np)
    {
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
            if (np > 0) {
                auto& aos = particle_tile.GetArrayOfStructs();
                dp = reinterpret_cast<ParticleReal*>(aos.data());
            } else {
                dp = nullptr;
            }
        } else {
            np = 0;
            dp = nullptr;
        }
    }

    void pt_amrex_fi_add_particle_i(FParticleContainer* particlecontainer,
                                 int lev, int grid, int tile, FParticleContainer::ParticleType* p)
    {
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto& particle_tile  = particle_level[std::make_pair(grid, tile)];
        particle_tile.push_back(*p);
    }

    void pt_amrex_fi_num_particles_i(FParticleContainer* particlecontainer,
                                  int lev, int grid, int tile, Long& np)
    {
        auto& particle_level = particlecontainer->GetParticles(lev);
        auto search = particle_level.find(std::make_pair(grid, tile));
        if (search != particle_level.end()) {
            auto& particle_tile = search->second;
            np = particle_tile.numParticles();
        } else {
            np = 0;
        }
    }
}
