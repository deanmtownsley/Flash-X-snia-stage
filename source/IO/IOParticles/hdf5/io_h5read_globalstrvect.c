/*This file allows us to get particle names to determine how many properties we
  have*/

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


void FTOC(io_h5read_globalstrvect)(hid_t* file_identifier,
                                   char strvect[][OUTPUT_PROP_LENGTH],
                                   char dataset_name[],
                                    int* num_elts){

  hid_t dataset,dataspace;
  herr_t status;
  
  hsize_t dimens_2d[2], maxdimens_2d[2];
  int i,  rank;
  hsize_t start_2d[2], count_2d[2], stride_2d[2];
  
  hid_t string_type;
  int string_size;
  
  /*manually set the string size */
  string_size = OUTPUT_PROP_LENGTH;
  
  /*setup the datatype for this string length */
  string_type = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(string_type, string_size);
  if(status < 0){
    Driver_abortC("Error: string\n");
  }

#ifdef FLASH_IO_ASYNC_HDF5
  dataset = H5Dopen_async(*file_identifier, dataset_name, H5P_DEFAULT);
#else
  dataset = H5Dopen(*file_identifier, dataset_name);
#endif
  dataspace = H5Dget_space(dataset);
  
  H5Sget_simple_extent_dims(dataspace, dimens_2d, maxdimens_2d);
  
  printf("Trying to read integer vector in dataset \"%s\", length in file is %d, number of elements requested %d\n", dataset_name, dimens_2d[0], *num_elts);

  if (dimens_2d[0] > *num_elts) {
    printf("Error: Integer vector in dataset \"%s\" in the file appears too long: %d > %d\n", dataset_name, dimens_2d[0], *num_elts);
    Driver_abortC("Error: Integer vector in the file appears too long\n");
  }
  *num_elts = dimens_2d[0];

  rank = 2;

  start_2d[0] = 0;
  start_2d[1] = 0;
  
  stride_2d[0] = 1;
  stride_2d[1] = 1;
  
  count_2d[0] = num_elts;
  count_2d[1] = 0;
  
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d, 
                        stride_2d, count_2d, NULL);
  
  if (status < 0){
    printf("Error: Unable to select hyperslab for global string vector \"%s\"\n", dataset_name);
    Driver_abortC("Error: Unable to select hyperslab for particle names\n");
  }
  
#ifdef FLASH_IO_ASYNC_HDF5
  status = H5Dread_async(dataset, string_type, H5S_ALL, H5S_ALL, 
                H5P_DEFAULT, strvect, io_es_id);
#else
  status = H5Dread(dataset, string_type, H5S_ALL, H5S_ALL, 
                H5P_DEFAULT, strvect);
#endif


  if (status < 0){
    printf("Error: Unable to read dataset for global string vector %d\n", status);
    Driver_abortC("Error: Unable to read dataset for global string vector\n");
  }

  H5Sclose(dataspace);
#ifdef FLASH_IO_ASYNC_HDF5
   H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif    
  
  return;
}
