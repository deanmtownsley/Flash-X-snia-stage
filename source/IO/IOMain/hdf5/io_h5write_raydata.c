#include <stdio.h>
#include "mangle_names.h"
#include "constants.h"
#include <hdf5.h>
#include <stdlib.h>
#include "hdf5_flash.h"
#include "Simulation.h"

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

void FTOC(io_h5write_raydata)(hid_t * pFileID,
			      int * totalPos,
			      int *ierr,
			      int *startPos,
			      int *npos,
			      double tags[],
			      double xpos[],
			      double ypos[], 
			      double zpos[], 
			      double power[],
			      int *rank)
{
  /* Open the ray dataset, it should already exist */
#ifdef FLASH_IO_ASYNC_HDF5
  hid_t dset = H5Dopen_async(*pFileID, "RayData", H5P_DEFAULT, io_es_id);
#else
  hid_t dset = H5Dopen(*pFileID, "RayData", H5P_DEFAULT);
#endif  
  if(dset < 0) {
    *ierr = dset;
    return;
  }

  /* Find out how big the dataset already is... */
  hsize_t dims[2], ndims;
  hid_t filespace = H5Dget_space (dset);
  if(filespace < 0) {
    *ierr = filespace;
    return;
  }

  *ierr = H5Sget_simple_extent_dims (filespace, dims, NULL);
  if(*ierr < 0) {
    return;
  }

  H5Sclose(filespace);

  /* Extend the dataset, add totalPos more entries... */
  hsize_t new_dims[2] = {dims[0] + *totalPos, 5};
  *ierr = H5Dextend (dset, new_dims);

  /* Dataset has been extended. Now, the data has to be written */
  filespace = H5Dget_space(dset);
  hsize_t offset[2] = {dims[0] + *startPos, 0};
  if(filespace < 0) {
    *ierr = filespace;
    return;
  }
  
  hsize_t slab_size[2] = {*npos, 5};

  /* printf("%i START_POS        = %i\n", *rank, *startPos); */
  /* printf("%i OLD DATASET SIZE = %i\t%i\n", *rank, dims[0], dims[1]); */
  /* printf("%i NEW DATASET SIZE = %i\t%i\n", *rank, new_dims[0], new_dims[1]); */
  /* printf("%i HYPERSLAB OFFSET = %i\t%i\n", *rank, offset[0], offset[1]); */
  /* printf("%i HYPERSLAB SIZE   = %i\t%i\n\n", *rank, slab_size[0], slab_size[1]); */

  *ierr = H5Sselect_hyperslab (filespace, H5S_SELECT_SET, offset, NULL,
			       slab_size, NULL);
  if(*ierr < 0) return;

  /* Define memory space */
  hsize_t max_size[2] = {H5S_UNLIMITED, 5};
  hid_t dataspace = H5Screate_simple (2, slab_size, max_size);
  if(dataspace < 0) {
    *ierr = dataspace;
    return;
  }

  /* Write the data to the hyperslab */
  double *data = (double*) malloc( (*npos)*5 * sizeof(double));

  int i;
  for(i = 0; i < *npos; i++) {
    data[i*5 + 0] = tags[i];
    data[i*5 + 1] = xpos[i];
    data[i*5 + 2] = ypos[i];
    data[i*5 + 3] = zpos[i];
    data[i*5 + 4] = power[i];
  }

  hid_t plist_id = H5Pcreate (H5P_DATASET_XFER);
#ifdef H5_HAVE_PARALLEL
    if (HDF5_MODE == COLLECTIVE) {
      *ierr = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
      if(*ierr < 0) return;
    } else {
      *ierr = H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_INDEPENDENT);
      if(*ierr < 0) return;
    }
#endif

#ifdef FLASH_IO_ASYNC_HDF5
  *ierr = H5Dwrite_async (dset, H5T_NATIVE_DOUBLE, dataspace, filespace,
          plist_id, data, io_es_id);
#else
  *ierr = H5Dwrite (dset, H5T_NATIVE_DOUBLE, dataspace, filespace,
  		    plist_id, data);
#endif

  free(data);
				  

  H5Pclose(plist_id);
  H5Sclose(dataspace);
  H5Sclose(filespace);
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async(dset, io_es_id);
#else
  H5Dclose(dset);
#endif  
}
