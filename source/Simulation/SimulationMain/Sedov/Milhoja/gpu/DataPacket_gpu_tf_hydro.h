#ifndef DATAPACKET_GPU_TF_HYDRO_UNIQUE_IFNDEF_H_
#define DATAPACKET_GPU_TF_HYDRO_UNIQUE_IFNDEF_H_

#include <Milhoja.h>
#include <Milhoja_real.h>
#include <Milhoja_DataPacket.h>
#include <climits>

using real = milhoja::Real;
using milhoja::FArray4D;
using milhoja::IntVect;
using milhoja::RealVect;

class DataPacket_gpu_tf_hydro : public milhoja::DataPacket {
public:
    // constructor
    DataPacket_gpu_tf_hydro(
    real dt
    
    );
    // destructor
    ~DataPacket_gpu_tf_hydro(void);

    //helper methods from base DataPacket class.
    std::unique_ptr<milhoja::DataPacket> clone(void) const override;
    DataPacket_gpu_tf_hydro(DataPacket_gpu_tf_hydro&) = delete;
    DataPacket_gpu_tf_hydro(const DataPacket_gpu_tf_hydro&) = delete;
    DataPacket_gpu_tf_hydro(DataPacket_gpu_tf_hydro&& packet) = delete;
    DataPacket_gpu_tf_hydro& operator=(DataPacket_gpu_tf_hydro&) = delete;
    DataPacket_gpu_tf_hydro& operator=(const DataPacket_gpu_tf_hydro&) = delete;
    DataPacket_gpu_tf_hydro& operator=(DataPacket_gpu_tf_hydro&& rhs) = delete;

    // pack and unpack functions from base class.
    void pack(void) override;
    void unpack(void) override;

    // TODO: Streams should be stored inside of an array.
    int extraAsynchronousQueue(const unsigned int id) override;
    void releaseExtraQueue(const unsigned int id) override;
    

    // DataPacket members are made public so a matching task function can easily access them.
    // Since both files are auto-generated and not maintained by humans, this is fine.
    real _dt_h;
    real* _dt_d;
    int _nTiles_h;
    int* _nTiles_d;
    real* _tile_deltas_d;
    int* _tile_lo_d;
    int* _tile_hi_d;
    int* _tile_lbound_d;
    real* _U_d;
    real* _U_p;
    real* _hydro_op1_auxc_d;
    real* _hydro_op1_flX_d;
    real* _hydro_op1_flY_d;
    real* _hydro_op1_flZ_d;
    
private:
    static constexpr std::size_t ALIGN_SIZE=16;
    static constexpr std::size_t pad(const std::size_t size) {
        return (((size + ALIGN_SIZE - 1) / ALIGN_SIZE) * ALIGN_SIZE);
    }

    // TODO: Streams should be stored inside of an array. Doing so would simplify the code
    // generation & source code for the stream functions.
    milhoja::Stream stream2_;
    milhoja::Stream stream3_;
    

    static constexpr std::size_t SIZE_DT = sizeof(real);
    static constexpr std::size_t SIZE_NTILES = sizeof(int);
    static constexpr std::size_t SIZE_TILE_DELTAS = MILHOJA_MDIM * sizeof(real);
    static constexpr std::size_t SIZE_TILE_LO = MILHOJA_MDIM * sizeof(int);
    static constexpr std::size_t SIZE_TILE_HI = MILHOJA_MDIM * sizeof(int);
    static constexpr std::size_t SIZE_TILE_LBOUND = MILHOJA_MDIM * sizeof(int);
    static constexpr std::size_t SIZE_U = (16 + 2 * 1 * MILHOJA_K1D) * (16 + 2 * 1 * MILHOJA_K2D) * (16 + 2 * 1 * MILHOJA_K3D) * (9 + 1 - 0) * sizeof(real);
    static constexpr std::size_t SIZE_HYDRO_OP1_AUXC = (18) * (18) * (18) * (1) * sizeof(real);
    static constexpr std::size_t SIZE_HYDRO_OP1_FLX = (19) * (18) * (18) * (5) * sizeof(real);
    static constexpr std::size_t SIZE_HYDRO_OP1_FLY = (18) * (19) * (18) * (5) * sizeof(real);
    static constexpr std::size_t SIZE_HYDRO_OP1_FLZ = (18) * (18) * (19) * (5) * sizeof(real);
    
};

#endif
