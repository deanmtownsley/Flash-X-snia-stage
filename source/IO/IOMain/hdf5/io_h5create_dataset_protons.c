#include <stdio.h>
#include "mangle_names.h"
#include "constants.h"
#include <hdf5.h>

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

void FTOC(io_h5create_dataset_protons)(const hid_t * const pFileID)
{
  /* Create the data space with unlimited 1st dimension,
     since we don't know how many protons will be written
     to the plot file. The 2nd dimension contains the
     proton tags and the x,y,z position coordinates. */

  hsize_t dims[2]       = {0, 4};
  hsize_t maxdims[2]    = {H5S_UNLIMITED, 4};
  hid_t   dataspaceID   = H5Screate_simple (2, dims, maxdims);
  
  /* Modify dataset creation properties, i.e. enable chunking.
     For this get a property identifier handle and use it to
     set chunking on it. Currently chunks of 256 are used. */

  hsize_t chunk_dims[2] = {256, 4};
  hid_t   propertyID    = H5Pcreate (H5P_DATASET_CREATE);
  herr_t  chunkError    = H5Pset_chunk (propertyID, 2, chunk_dims);
  
  /* Create a new dataset within the file using the dataspace
     and property identifiers. */
#ifdef FLASH_IO_ASYNC_HDF5
  hid_t   datasetID     = H5Dcreate_async (*pFileID,             /* file identifier pointer     */
                                     "ProtonData",         /* name of dataset created     */
                                     H5T_NATIVE_DOUBLE,    /* data type handle            */
                                     dataspaceID,          /* dataspace identifier handle */
			             H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT, io_es_id);          /* property identifier handle  */

#else
  hid_t   datasetID     = H5Dcreate (*pFileID,             /* file identifier pointer     */
                                     "ProtonData",         /* name of dataset created     */
                                     H5T_NATIVE_DOUBLE,    /* data type handle            */
                                     dataspaceID,          /* dataspace identifier handle */
			             H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);          /* property identifier handle  */

#endif
  /* Return the identifier handles. */

  H5Sclose (dataspaceID);
  H5Pclose (propertyID);
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async (datasetID, io_es_id);
#else
  H5Dclose (datasetID);
#endif  
}
