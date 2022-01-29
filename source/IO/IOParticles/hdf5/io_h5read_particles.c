#include "mangle_names.h"
#include <assert.h>
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "Simulation.h"
#include "constants.h"

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

int Driver_abortC(char* message);



/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(get_numparticles)(hid_t *file_identifier,
                   int *numParticles)
{
  
  /* 
   * get_numParticles: look at the particles entry and get the size
   *                of the dataset based on the dimensions of the array 
   *
   * input:  file handle          (integer)
   *
   * output: number of particles  (integer)
   *
   * return value: -1 indicates an error getting the dimension
   *
   */

  herr_t   status;

  hsize_t maximum_dims[10];
  hsize_t dataspace_dims[10];

  hid_t dataset, dataspace;

  herr_t (*old_func)(void*);
  void *old_client_data;
   
  /* temporarily turn off error handling, to probe for the file's existence */
#ifdef FLASH_IO_ASYNC_HDF5
  H5Eget_auto1(&old_func, &old_client_data);
  H5Eset_auto1(NULL, NULL);
#else
  H5Eget_auto(H5P_DEFAULT, &old_func, &old_client_data);
  H5Eset_auto(H5P_DEFAULT, NULL, NULL);
#endif
  /* get the dataset ID for the coordinates record, and the dataspace in
     this record */
#ifdef FLASH_IO_ASYNC_HDF5
  dataset = H5Dopen_async(*file_identifier, "particle tracers", H5P_DEFAULT, io_es_id);
#else
  dataset = H5Dopen(*file_identifier, "particle tracers", H5P_DEFAULT);
#endif     

  /* restore the error handling */
#ifdef FLASH_IO_ASYNC_HDF5
  H5Eset_auto1(old_func, old_client_data);
#else
  H5Eset_auto(H5P_DEFAULT, old_func, old_client_data);
#endif  

  if (dataset >= 0) {
    dataspace = H5Dget_space(dataset);

    /* get the dimensions of the dataspace */
    status = H5Sget_simple_extent_dims(dataspace, dataspace_dims, maximum_dims);
    if(status < 0){
      printf("No particles found in particle tracers dataset!\n");
      Driver_abortC("No particles found in particle tracers dataset!\n");
    }

    *numParticles = dataspace_dims[0];

#ifdef FLASH_IO_ASYNC_HDF5
    H5Dclose_async(dataset, io_es_id);  
#else
    H5Dclose(dataset);  
#endif
    H5Sclose(dataspace);

  } else {
    *numParticles = 0;
  }

}

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5read_particles)(hid_t* file_identifier,
                         double* particles,
                         int* localnp,
                         int* npart_props,
                         int* particle_offset)
{

  hid_t    dataspace, dataset, memspace, xfer_plist;
  herr_t   status, herr;

  int      rank;
  hsize_t  dimens_2d[2], start_2d[2], count_2d[2], stride_2d[2];


  xfer_plist = H5Pcreate(H5P_DATASET_XFER);
#ifdef IO_HDF5_PARALLEL
  if (HDF5_MODE == COLLECTIVE) {
    herr = H5Pset_dxpl_mpio(xfer_plist, H5FD_MPIO_COLLECTIVE);
    assert(herr >= 0);
  }
#endif


  /* open the dataset */
#ifdef FLASH_IO_ASYNC_HDF5
  dataset = H5Dopen_async(*file_identifier, "tracer particles", H5P_DEFAULT, io_es_id);
#else
  dataset = H5Dopen(*file_identifier, "tracer particles", H5P_DEFAULT);
#endif  
  if(dataset < 0){
    Driver_abortC("tracer particles not found\n");
  }

  /* get the particle data space */
  dataspace  = H5Dget_space(dataset);


  /* select this processor's hyperslab */
  start_2d[0] = (hsize_t) (*particle_offset);
  start_2d[1] = (hsize_t) 0;

  stride_2d[0] = 1;
  stride_2d[1] = 1;

  count_2d[0]  = (hsize_t) (*localnp);
  count_2d[1]  = (hsize_t) (*npart_props);

  status     = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_2d,
                                   stride_2d, count_2d, NULL);
  if (status < 0){
    printf("Error: Unable to select hyperslab for particles dataspace\n");
    Driver_abortC("Error: Unable to select hyperslab for particles dataspace\n");
  }

  if(*localnp > 0){
    
    rank       = 2;
    dimens_2d[0] = (*localnp);
    dimens_2d[1] = (*npart_props);
    
    memspace   = H5Screate_simple(rank, dimens_2d, NULL);
    
    /* read the data */
#ifdef FLASH_IO_ASYNC_HDF5
    status = H5Dread_async(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                 xfer_plist, particles, io_es_id);
#else
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                 xfer_plist, particles);
#endif    

    
    if (status < 0){
      printf("Error: Unable to read particles from dataset\n");
      Driver_abortC("Error: Unable to read particles from dataset\n");
    }


  /* release resources */
    status = H5Sclose(memspace);

  }

  H5Pclose(xfer_plist);
  status = H5Sclose(dataspace);
#ifdef FLASH_IO_ASYNC_HDF5
   H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif
}

