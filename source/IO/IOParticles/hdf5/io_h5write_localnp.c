#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "Simulation.h"
#include "constants.h"

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif


int Driver_abortC(char* message);

/* xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx */

void FTOC(io_h5write_localnp)(int* myPE, 
                        hid_t* file_identifier,
                        int* localnp,
                        int* local_blocks,
                        int* total_blocks,
                        int* global_offset)
{
  hid_t dataspace, dataset, memspace, dataset_plist, dxfer_template;
  herr_t status;

  int rank;
  hsize_t dimens_1d;

  hsize_t start_1d;
  hsize_t stride_1d, count_1d;
  /* set the dimensions of the dataset */
  rank = 1;
  dimens_1d = *total_blocks;
 
  start_1d = (hsize_t) (*global_offset);
  stride_1d = 1;
  count_1d = (hsize_t) (*local_blocks);

 
  dataspace = H5Screate_simple(rank, &dimens_1d, NULL);
  if(dataspace < 0) {
     Driver_abortC("Error: H5Screate_simple io_h5write_localnp\n");
  }


  /* create the hyperslab -- this will differ on the different processors */
  status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, &start_1d, 
                        &stride_1d, &count_1d, NULL);
  if(status < 0) {
     Driver_abortC("Error: H5Sselect_hyperslab io_h5write_localnp\n");
  }

 
  /* create the memory space */
  dimens_1d = *local_blocks;
  memspace = H5Screate_simple(rank, &dimens_1d, NULL);
  if(memspace < 0) {
     Driver_abortC("Error: H5Screate_simple mem io_h5write_localnp\n");
  }
  


  /*DEV: this line was only in the serial version */
  dataset_plist = H5Pcreate(H5P_DATASET_CREATE);
  if(dataset_plist < 0) {
     Driver_abortC("Error: dataset_plist io_h5write_localnp\n");
  }


  /*This part is necessary when serial IO is used
    The master proc writes all the data but only creates
    the dataset once.  In the parallel IO case each 
    proc calls H5Dcreate*/
  if ((*myPE == MASTER_PE) && (*global_offset != 0)) {
#ifdef FLASH_IO_ASYNC_HDF5
     dataset = H5Dopen_async(*file_identifier, "localnp", H5P_DEFAULT, io_es_id); 
#else
     dataset = H5Dopen(*file_identifier, "localnp", H5P_DEFAULT); 
#endif     
    if(dataset < 0) {
       Driver_abortC("Error: H5Dopen io_h5write_localnp\n");
    }
   
  }else {

     /* create the dataset */
#ifdef FLASH_IO_ASYNC_HDF5
    dataset = H5Dcreate_async(*file_identifier, "localnp", H5T_NATIVE_INT,
                  dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT, io_es_id);
#else
    dataset = H5Dcreate(*file_identifier, "localnp", H5T_NATIVE_INT,
                  dataspace, H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
#endif     

    if(dataset < 0) {
       Driver_abortC("Error: H5Dcreate io_h5write_localnp\n");
    }

    /*dataset_plist was H5P_DEFAULT*/
  }  
  dxfer_template = H5Pcreate(H5P_DATASET_XFER);
#ifdef IO_HDF5_PARALLEL
  if(HDF5_MODE == COLLECTIVE)
    H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
#endif
  
  /* write the data */
#ifdef FLASH_IO_ASYNC_HDF5
  status = H5Dwrite_async(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                 dxfer_template, localnp, io_es_id);
#else
  status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, 
                dxfer_template, localnp);
#endif  

  if(status < 0) {
     Driver_abortC("Error: H5Dwrite io_h5write_localnp\n");
  }


  H5Pclose(dxfer_template);
  H5Pclose(dataset_plist);
  H5Sclose(memspace); 
  H5Sclose(dataspace);
  
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif  
  
}























