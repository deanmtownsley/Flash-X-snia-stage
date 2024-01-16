/* _connector:get_host_members */
const int queue1_h = packet_h->asynchronousQueue();
const int _nTiles_h = packet_h->_nTiles_h;

/* _connector:c2f_argument_list */
void* packet_h,
const int queue1_h,
const int _nTiles_h,
const void* _nTiles_d,
const void* _external_hydro_op1_dt_d,
const void* _tile_deltas_d,
const void* _tile_lo_d,
const void* _tile_hi_d,
const void* _tile_lbound_d,
const void* _CC_1_d,
const void* _scratch_hydro_op1_auxC_d,
const void* _scratch_hydro_op1_flX_d,
const void* _scratch_hydro_op1_flY_d,
const void* _scratch_hydro_op1_flZ_d

/* _connector:c2f_arguments */
packet_h,
queue1_h,
_nTiles_h,
_nTiles_d,
_external_hydro_op1_dt_d,
_tile_deltas_d,
_tile_lo_d,
_tile_hi_d,
_tile_lbound_d,
_CC_1_d,
_scratch_hydro_op1_auxC_d,
_scratch_hydro_op1_flX_d,
_scratch_hydro_op1_flY_d,
_scratch_hydro_op1_flZ_d

/* _connector:get_device_members */
void* _nTiles_d = static_cast<void*>( packet_h->_nTiles_d );
void* _external_hydro_op1_dt_d = static_cast<void*>( packet_h->_external_hydro_op1_dt_d );
void* _tile_deltas_d = static_cast<void*>( packet_h->_tile_deltas_d );
void* _tile_lo_d = static_cast<void*>( packet_h->_tile_lo_d );
void* _tile_hi_d = static_cast<void*>( packet_h->_tile_hi_d );
void* _tile_lbound_d = static_cast<void*>( packet_h->_tile_lbound_d );
void* _CC_1_d = static_cast<void*>( packet_h->_CC_1_d );
void* _scratch_hydro_op1_auxC_d = static_cast<void*>( packet_h->_scratch_hydro_op1_auxC_d );
void* _scratch_hydro_op1_flX_d = static_cast<void*>( packet_h->_scratch_hydro_op1_flX_d );
void* _scratch_hydro_op1_flY_d = static_cast<void*>( packet_h->_scratch_hydro_op1_flY_d );
void* _scratch_hydro_op1_flZ_d = static_cast<void*>( packet_h->_scratch_hydro_op1_flZ_d );

/* _connector:instance_args */
real external_hydro_op1_dt,void** packet

/* _connector:host_members */
external_hydro_op1_dt