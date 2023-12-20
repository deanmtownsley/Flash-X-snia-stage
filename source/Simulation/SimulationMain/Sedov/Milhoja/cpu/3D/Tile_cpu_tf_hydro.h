#ifndef MILHOJA_GENERATED_TILE_CPU_TF_HYDRO_H__
#define MILHOJA_GENERATED_TILE_CPU_TF_HYDRO_H__

#include <Milhoja_real.h>
#include <Milhoja_TileWrapper.h>

struct Tile_cpu_tf_hydro : public milhoja::TileWrapper {
    Tile_cpu_tf_hydro(const milhoja::Real external_hydro_op1_dt,
                      const int external_hydro_op1_eosMode);
    ~Tile_cpu_tf_hydro(void);

    Tile_cpu_tf_hydro(Tile_cpu_tf_hydro&)                  = delete;
    Tile_cpu_tf_hydro(const Tile_cpu_tf_hydro&)            = delete;
    Tile_cpu_tf_hydro(Tile_cpu_tf_hydro&&)                 = delete;
    Tile_cpu_tf_hydro& operator=(Tile_cpu_tf_hydro&)       = delete;
    Tile_cpu_tf_hydro& operator=(const Tile_cpu_tf_hydro&) = delete;
    Tile_cpu_tf_hydro& operator=(Tile_cpu_tf_hydro&&)      = delete;

    std::unique_ptr<milhoja::TileWrapper> clone(std::shared_ptr<milhoja::Tile>&& tileToWrap) const override;

    milhoja::Real  external_hydro_op1_dt_;
    int            external_hydro_op1_eosMode_;

    static void acquireScratch(void);
    static void releaseScratch(void);

    constexpr static std::size_t  SCRATCH_HYDRO_OP1_AUXC_SIZE_ =
                      18
                    * 18
                    * 18;
    constexpr static std::size_t  SCRATCH_HYDRO_OP1_FLX_SIZE_ =
                      17
                    * 16
                    * 16
                    *  5;
    constexpr static std::size_t  SCRATCH_HYDRO_OP1_FLY_SIZE_ =
                      16
                    * 17
                    * 16
                    *  5;
    constexpr static std::size_t  SCRATCH_HYDRO_OP1_FLZ_SIZE_ =
                      16
                    * 16
                    * 17
                    *  5;

    static void* scratch_hydro_op1_auxc_;
    static void* scratch_hydro_op1_flx_;
    static void* scratch_hydro_op1_fly_;
    static void* scratch_hydro_op1_flz_;
};

#endif
