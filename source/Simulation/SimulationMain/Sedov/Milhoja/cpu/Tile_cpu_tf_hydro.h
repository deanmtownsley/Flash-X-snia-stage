#ifndef MILHOJA_GENERATED_TILE_CPU_TF_HYDRO_H__
#define MILHOJA_GENERATED_TILE_CPU_TF_HYDRO_H__

#include <Milhoja_TileWrapper.h>

struct Tile_cpu_tf_hydro : public milhoja::TileWrapper {
    Tile_cpu_tf_hydro(const milhoja::Real dt);
    ~Tile_cpu_tf_hydro(void);

    Tile_cpu_tf_hydro(Tile_cpu_tf_hydro&)                  = delete;
    Tile_cpu_tf_hydro(const Tile_cpu_tf_hydro&)            = delete;
    Tile_cpu_tf_hydro(Tile_cpu_tf_hydro&&)                 = delete;
    Tile_cpu_tf_hydro& operator=(Tile_cpu_tf_hydro&)       = delete;
    Tile_cpu_tf_hydro& operator=(const Tile_cpu_tf_hydro&) = delete;
    Tile_cpu_tf_hydro& operator=(Tile_cpu_tf_hydro&&)      = delete;

    std::unique_ptr<milhoja::TileWrapper> clone(std::shared_ptr<milhoja::Tile>&& tileToWrap) const override;

    milhoja::Real  dt_;

    static void acquireScratch(void);
    static void releaseScratch(void);

    // TODO: This is large enough that it will work for 2D and 3D.  This should
    // be dimension-specific once we have actual on-demand code generation.
    constexpr static std::size_t  HYDRO_OP1_AUXC_SIZE_ =
                      18
                    * 18
                    * 18;

    static void* hydro_op1_auxc_;
};

#endif
