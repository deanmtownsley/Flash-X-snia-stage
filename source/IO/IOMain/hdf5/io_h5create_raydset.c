#include <stdio.h>
#include "mangle_names.h"
#include "constants.h"
#include <hdf5.h>
#include "Simulation.h"

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

void FTOC(io_h5create_raydset)(const hid_t * const pFileID)
{
  /* Create the data space with unlimited dimensions. */
  hsize_t dims[2] = {0, 5};
  hsize_t maxdims[2] = {H5S_UNLIMITED, 5};
  hid_t dataspace = H5Screate_simple (2, dims, maxdims);
  
  /* Modify dataset creation properties, i.e. enable chunking  */
  hsize_t chunk_dims[2] = {256, 5}; /* TODO: Pass in good chunk size estimate */
  hid_t cparms = H5Pcreate (H5P_DATASET_CREATE);
  H5Pset_chunk ( cparms, 2, chunk_dims);
  
  /* Create a new dataset within the file using cparms
     creation properties.  */
#ifdef FLASH_IO_ASYNC_HDF5
  hid_t dsetid = H5Dcreate_async (*pFileID, "RayData", H5T_NATIVE_DOUBLE, dataspace,
			    H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT, io_es_id);
#else
  hid_t dsetid = H5Dcreate (*pFileID, "RayData", H5T_NATIVE_DOUBLE, dataspace,
			    H5P_DEFAULT,H5P_DEFAULT, H5P_DEFAULT);
#endif


  H5Sclose(dataspace);
  H5Pclose(cparms);
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async(dsetid, io_es_id);
#else
  H5Dclose(dsetid);
#endif

}
