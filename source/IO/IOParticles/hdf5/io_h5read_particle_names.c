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


void FTOC(io_h5read_particle_names)(hid_t* file_identifier,
                                   char part_names[][OUTPUT_PROP_LENGTH],
                                    int* num_part_props){

  hid_t dataset;
  herr_t status;
  
  hsize_t dimens_2d[2], maxdimens_2d[2];
  int i;
  
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
  dataset = H5Dopen_async(*file_identifier, "particle names", H5P_DEFAULT, io_es_id);
#else
  dataset = H5Dopen(*file_identifier, "particle names", H5P_DEFAULT);
#endif
  if (status < 0){
    printf("Error: Unable to open dataset for particle names %d\n", status);
    Driver_abortC("Error: Unable to open dataset for particle names\n");
  }

  
  status = H5Dread(dataset, string_type, H5S_ALL, H5S_ALL, 
                H5P_DEFAULT, part_names);

  if (status < 0){
    printf("Error: Unable to read dataset for particle names %d\n", status);
    Driver_abortC("Error: Unable to read dataset for particle names\n");
  }

#ifdef FLASH_IO_ASYNC_HDF5
   H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif

  if(status < 0){
    Driver_abortC("Error: Close dataset\n");
  }

  status = H5Tclose(string_type);
  if(status < 0){
    Driver_abortC("Error: Free string type\n");
  }

  return;
}
