/* _connector:size_constructor */
SIZE_EXTERNAL_HYDRO_OP1_DT
 + SIZE_NTILES

/* _connector:size_tile_metadata */
SIZE_TILE_DELTAS
 + SIZE_TILE_LO
 + SIZE_TILE_HI
 + SIZE_TILE_LBOUND

/* _connector:size_tile_in */
0

/* _connector:size_tile_in_out */
SIZE_CC_1

/* _connector:size_tile_out */
0

/* _connector:size_tile_scratch */
SIZE_SCRATCH_HYDRO_OP1_AUXC
 + SIZE_SCRATCH_HYDRO_OP1_FLX
 + SIZE_SCRATCH_HYDRO_OP1_FLY
 + SIZE_SCRATCH_HYDRO_OP1_FLZ

/* _connector:constructor_args */
real external_hydro_op1_dt

/* _connector:set_members */
_external_hydro_op1_dt_h{external_hydro_op1_dt},
_external_hydro_op1_dt_d{nullptr},
_nTiles_h{0},
_nTiles_d{nullptr},
_tile_deltas_d{nullptr},
_tile_lo_d{nullptr},
_tile_hi_d{nullptr},
_tile_lbound_d{nullptr},
_CC_1_d{nullptr},
_CC_1_p{nullptr},
_scratch_hydro_op1_auxC_d{nullptr},
_scratch_hydro_op1_flX_d{nullptr},
_scratch_hydro_op1_flY_d{nullptr},
_scratch_hydro_op1_flZ_d{nullptr}

/* _connector:host_members */
_external_hydro_op1_dt_h

/* _connector:size_determination */
static constexpr std::size_t SIZE_EXTERNAL_HYDRO_OP1_DT = sizeof(real);
static constexpr std::size_t SIZE_NTILES = sizeof(int);
static constexpr std::size_t SIZE_TILE_DELTAS = MILHOJA_MDIM * sizeof(real);
static constexpr std::size_t SIZE_TILE_LO = MILHOJA_MDIM * sizeof(int);
static constexpr std::size_t SIZE_TILE_HI = MILHOJA_MDIM * sizeof(int);
static constexpr std::size_t SIZE_TILE_LBOUND = MILHOJA_MDIM * sizeof(int);
static constexpr std::size_t SIZE_CC_1 = (8 + 2 * 1 * MILHOJA_K1D) * (8 + 2 * 1 * MILHOJA_K2D) * (1 + 2 * 1 * MILHOJA_K3D) * (8 + 1 - 0) * sizeof(real);
static constexpr std::size_t SIZE_SCRATCH_HYDRO_OP1_AUXC = (10) * (10) * (1) * sizeof(real);
static constexpr std::size_t SIZE_SCRATCH_HYDRO_OP1_FLX = (11) * (10) * (1) * (5) * sizeof(real);
static constexpr std::size_t SIZE_SCRATCH_HYDRO_OP1_FLY = (10) * (11) * (1) * (5) * sizeof(real);
static constexpr std::size_t SIZE_SCRATCH_HYDRO_OP1_FLZ = (1) * (1) * (1) * (1) * sizeof(real);

/* _connector:pointers_constructor */
real* _external_hydro_op1_dt_p = static_cast<real*>( static_cast<void*>(ptr_p) );
_external_hydro_op1_dt_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_p+=SIZE_EXTERNAL_HYDRO_OP1_DT;
ptr_d+=SIZE_EXTERNAL_HYDRO_OP1_DT;

int* _nTiles_p = static_cast<int*>( static_cast<void*>(ptr_p) );
_nTiles_d = static_cast<int*>( static_cast<void*>(ptr_d) );
ptr_p+=SIZE_NTILES;
ptr_d+=SIZE_NTILES;


/* _connector:memcpy_constructor */
std::memcpy(_external_hydro_op1_dt_p, static_cast<void*>(&_external_hydro_op1_dt_h), SIZE_EXTERNAL_HYDRO_OP1_DT);
std::memcpy(_nTiles_p, static_cast<void*>(&_nTiles_h), SIZE_NTILES);

/* _connector:nTiles_value */
// Check for overflow first to avoid UB
// TODO: Should casting be checked here or in base class?
if (tiles_.size() > INT_MAX)
	throw std::overflow_error("[_param:class_name pack] nTiles was too large for int.");
_nTiles_h = static_cast<int>(tiles_.size());
/* _connector:public_members */
real _external_hydro_op1_dt_h;
real* _external_hydro_op1_dt_d;
int _nTiles_h;
int* _nTiles_d;
real* _tile_deltas_d;
int* _tile_lo_d;
int* _tile_hi_d;
int* _tile_lbound_d;
real* _CC_1_d;
real* _CC_1_p;
real* _scratch_hydro_op1_auxC_d;
real* _scratch_hydro_op1_flX_d;
real* _scratch_hydro_op1_flY_d;
real* _scratch_hydro_op1_flZ_d;

/* _connector:pointers_tile_metadata */
real* _tile_deltas_p = static_cast<real*>( static_cast<void*>(ptr_p) );
_tile_deltas_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_p+=_nTiles_h * SIZE_TILE_DELTAS;
ptr_d+=_nTiles_h * SIZE_TILE_DELTAS;

int* _tile_lo_p = static_cast<int*>( static_cast<void*>(ptr_p) );
_tile_lo_d = static_cast<int*>( static_cast<void*>(ptr_d) );
ptr_p+=_nTiles_h * SIZE_TILE_LO;
ptr_d+=_nTiles_h * SIZE_TILE_LO;

int* _tile_hi_p = static_cast<int*>( static_cast<void*>(ptr_p) );
_tile_hi_d = static_cast<int*>( static_cast<void*>(ptr_d) );
ptr_p+=_nTiles_h * SIZE_TILE_HI;
ptr_d+=_nTiles_h * SIZE_TILE_HI;

int* _tile_lbound_p = static_cast<int*>( static_cast<void*>(ptr_p) );
_tile_lbound_d = static_cast<int*>( static_cast<void*>(ptr_d) );
ptr_p+=_nTiles_h * SIZE_TILE_LBOUND;
ptr_d+=_nTiles_h * SIZE_TILE_LBOUND;


/* _connector:memcpy_tile_metadata */
real _tile_deltas_h[MILHOJA_MDIM] = { deltas.I(), deltas.J(), deltas.K() };
char_ptr = static_cast<char*>(static_cast<void*>(_tile_deltas_p)) + n * SIZE_TILE_DELTAS;
std::memcpy(static_cast<void*>(char_ptr), static_cast<void*>(_tile_deltas_h), SIZE_TILE_DELTAS);

int _tile_lo_h[MILHOJA_MDIM] = { lo.I()+1, lo.J()+1, lo.K()+1 };
char_ptr = static_cast<char*>(static_cast<void*>(_tile_lo_p)) + n * SIZE_TILE_LO;
std::memcpy(static_cast<void*>(char_ptr), static_cast<void*>(_tile_lo_h), SIZE_TILE_LO);

int _tile_hi_h[MILHOJA_MDIM] = { hi.I()+1, hi.J()+1, hi.K()+1 };
char_ptr = static_cast<char*>(static_cast<void*>(_tile_hi_p)) + n * SIZE_TILE_HI;
std::memcpy(static_cast<void*>(char_ptr), static_cast<void*>(_tile_hi_h), SIZE_TILE_HI);

int _tile_lbound_h[MILHOJA_MDIM] = { loGC.I()+1, loGC.J()+1, loGC.K()+1 };
char_ptr = static_cast<char*>(static_cast<void*>(_tile_lbound_p)) + n * SIZE_TILE_LBOUND;
std::memcpy(static_cast<void*>(char_ptr), static_cast<void*>(_tile_lbound_h), SIZE_TILE_LBOUND);


/* _connector:tile_descriptor */
const auto deltas = tileDesc_h->deltas();
const auto lo = tileDesc_h->lo();
const auto hi = tileDesc_h->hi();
const auto loGC = tileDesc_h->loGC();

/* _connector:pointers_tile_in */

/* _connector:memcpy_tile_in */

/* _connector:pointers_tile_in_out */
_CC_1_p = static_cast<real*>( static_cast<void*>(ptr_p) );
_CC_1_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_p+=_nTiles_h * SIZE_CC_1;
ptr_d+=_nTiles_h * SIZE_CC_1;


/* _connector:memcpy_tile_in_out */
real* CC_1_d = tileDesc_h->dataPtr();
constexpr std::size_t offset_CC_1 = (8 + 2 * 1 * MILHOJA_K1D) * (8 + 2 * 1 * MILHOJA_K2D) * (1 + 2 * 1 * MILHOJA_K3D) * static_cast<std::size_t>(0);
constexpr std::size_t nBytes_CC_1 = (8 + 2 * 1 * MILHOJA_K1D) * (8 + 2 * 1 * MILHOJA_K2D) * (1 + 2 * 1 * MILHOJA_K3D) * ( 8 - 0 + 1 ) * sizeof(real);
char_ptr = static_cast<char*>( static_cast<void*>(_CC_1_p) ) + n * SIZE_CC_1;
std::memcpy(static_cast<void*>(char_ptr), static_cast<void*>(CC_1_d + offset_CC_1), nBytes_CC_1);


/* _connector:unpack_tile_in_out */
constexpr std::size_t offset_CC_1 = (8 + 2 * 1 * MILHOJA_K1D) * (8 + 2 * 1 * MILHOJA_K2D) * (1 + 2 * 1 * MILHOJA_K3D) * static_cast<std::size_t>(0);
real*        start_h_CC_1 = CC_1_data_h + offset_CC_1;
const real*  start_p_CC_1 = CC_1_data_p + offset_CC_1;
constexpr std::size_t nBytes_CC_1 = (8 + 2 * 1 * MILHOJA_K1D) * (8 + 2 * 1 * MILHOJA_K2D) * (1 + 2 * 1 * MILHOJA_K3D) * ( 7 - 0 + 1 ) * sizeof(real);
std::memcpy(static_cast<void*>(start_h_CC_1), static_cast<const void*>(start_p_CC_1), nBytes_CC_1);


/* _connector:in_pointers */
real* CC_1_data_h = tileDesc_h->dataPtr();

/* _connector:out_pointers */
real* CC_1_data_p = static_cast<real*>( static_cast<void*>( static_cast<char*>( static_cast<void*>( _CC_1_p ) ) + n * SIZE_CC_1 ) );

/* _connector:pointers_tile_out */

/* _connector:memcpy_tile_out */

/* _connector:unpack_tile_out */

/* _connector:pointers_tile_scratch */
_scratch_hydro_op1_auxC_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_d += _nTiles_h * SIZE_SCRATCH_HYDRO_OP1_AUXC;

_scratch_hydro_op1_flX_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_d += _nTiles_h * SIZE_SCRATCH_HYDRO_OP1_FLX;

_scratch_hydro_op1_flY_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_d += _nTiles_h * SIZE_SCRATCH_HYDRO_OP1_FLY;

_scratch_hydro_op1_flZ_d = static_cast<real*>( static_cast<void*>(ptr_d) );
ptr_d += _nTiles_h * SIZE_SCRATCH_HYDRO_OP1_FLZ;


/* _connector:memcpy_tile_scratch */

