#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
#include "Simulation.h"
#include "constants.h"

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

int Driver_abortC(char* message);


void FTOC(io_h5read_num_props)(hid_t* file_identifier,
                                    int* num_part_props){

  hid_t dataset,dataspace, memspace;
  herr_t status;

  hsize_t dimens_2d[2], maxdimens_2d[2];

#ifdef FLASH_IO_ASYNC_HDF5
  dataset = H5Dopen_async(*file_identifier, "particle names", H5P_DEFAULT, io_es_id);
#else
  dataset = H5Dopen(*file_identifier, "particle names", H5P_DEFAULT);
#endif
  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
  
  *num_part_props = dimens_2d[0];
  
  
  H5Sclose(dataspace);
#ifdef FLASH_IO_ASYNC_HDF5
   H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif  
  return;
}
