#include <stdio.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <assert.h>
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

/* This function writes out a single unknown (passed from the checkpoint 
   routine), giving the record a label from the varnames or species
   database */
void FTOC(io_h5read_unknowns)(hid_t* file_identifier,
                       int* nxb,             /* num of zones to store in x */
                   int* nyb,             /* num of zones to store in y */
                   int* nzb,             /* num of zones to store in z */
                   double* unknowns,     /* [mblk][NZB][NYB][NXB][var] */
                   char record_label[5], /* add char--null termination */
                   int* local_blocks,
                   int* total_blocks,
                   int* global_offset,
                   int* doread)
{
  hid_t dataspace, dataset, memspace, dxfer_template;
  herr_t status;

  int rank;
  hsize_t dimens_5d[5];

  hsize_t start_4d[4], start_5d[5];
  hsize_t stride_4d[4], stride_5d[5], count_4d[4], count_5d[5];

  char record_label_new[5];

  int ierr;
  herr_t herr;

  herr_t (*old_func)(void*) = NULL;
  void *old_client_data = NULL;
  size_t spaceTrimLen;

  dxfer_template = H5Pcreate(H5P_DATASET_XFER);
#ifdef IO_HDF5_PARALLEL
  if (HDF5_MODE == COLLECTIVE) {
    herr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_COLLECTIVE);
    assert(herr >= 0);
  }
#endif


  /* the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination */
  strncpy(record_label_new, record_label,4);

  *(record_label_new + 4) = '\0';


  /* open the dataset
     suppress the error message if the variable does not exist in file */
  herr = H5Eget_auto(H5P_DEFAULT, &old_func, &old_client_data);
  assert(herr >= 0);
  herr = H5Eset_auto(H5P_DEFAULT, NULL, NULL);
  assert(herr >= 0);
#ifdef FLASH_IO_ASYNC_HDF5
      dataset = H5Dopen_async(*file_identifier, record_label_new,H5P_DEFAULT, io_es_id); 
#else
      dataset = H5Dopen(*file_identifier, record_label_new,H5P_DEFAULT);
#endif  
  /* If the open fails then we attempt to open the dataset again with a
     space trimmed name.  We do this because type based I/O in FLASH
     creates datasets which have names that are stripped of space padding */
  if (dataset < 0) {
    spaceTrimLen = strcspn(record_label_new, " ");
    if (spaceTrimLen > 0) {
      *(record_label_new + spaceTrimLen) = '\0';

#ifdef FLASH_IO_ASYNC_HDF5
      dataset = H5Dopen_async(*file_identifier, record_label_new,H5P_DEFAULT, io_es_id); 
#else
      dataset = H5Dopen(*file_identifier, record_label_new,H5P_DEFAULT);
#endif
    }
  }
  herr = H5Eset_auto(H5P_DEFAULT, old_func, old_client_data);
  assert(herr >= 0);

  if (dataset < 0) {
    printf("couldn't find variable '%s' in the file, so skipping it\n", record_label_new);
  } else {
    dataspace = H5Dget_space(dataset);
    if(dataspace < 0) Driver_abortC("Error: negative return from dataspace H5Dget_space");

    if (*local_blocks > 0) {
      /* create the hyperslab -- this will differ on the different processors */
      start_4d[0] = (hsize_t) (*global_offset);
      start_4d[1] = 0;
      start_4d[2] = 0;
      start_4d[3] = 0;

      stride_4d[0] = 1;
      stride_4d[1] = 1;
      stride_4d[2] = 1;
      stride_4d[3] = 1;

      count_4d[0] = (hsize_t) (*local_blocks);
      count_4d[1] = *nzb;
      count_4d[2] = *nyb;
      count_4d[3] = *nxb;

      status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d, 
				   stride_4d, count_4d, NULL);
  
      if (status < 0){
	printf("Error: Unable to select hyperslab for unknowns dataspace\n");
	Driver_abortC("Error: Unable to select hyperslab for unknowns dataspace\n");
      }

      /* create the memory space */
      rank = 5;
      dimens_5d[0] = *local_blocks;
      dimens_5d[1] = *nzb;
      dimens_5d[2] = *nyb;
      dimens_5d[3] = *nxb;
      dimens_5d[4] = 1;
    } else {
      /* MPI rank owns zero blocks - create a memory space of
	 size 1 so that H5Screate_simple does not abort */
      rank = 5;
      dimens_5d[0] = 1;
      dimens_5d[1] = 1;
      dimens_5d[2] = 1;
      dimens_5d[3] = 1;
      dimens_5d[4] = 1;
    }
    memspace = H5Screate_simple(rank, dimens_5d, NULL);
    if(memspace < 0) Driver_abortC("Error: negative return from memspace H5Screate_simple");

    /*
    start_5d[0] = 0;
    start_5d[1] = (*nguard)*k3d;
    start_5d[2] = (*nguard)*k2d;
    start_5d[3] = (*nguard);
    start_5d[4] = *index-1;  

    stride_5d[0] = 1;
    stride_5d[1] = 1;
    stride_5d[2] = 1;
    stride_5d[3] = 1;
    stride_5d[4] = 1;

    count_5d[0] = *local_blocks;
    count_5d[1] = *nzb;
    count_5d[2] = *nyb;
    count_5d[3] = *nxb;
    count_5d[4] = 1;
   
    ierr = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
                         start_5d, stride_5d, count_5d, NULL);
    

    if (ierr < 0){
      printf("Error: Unable to select hyperslab for coordinates memspace\n");
      Driver_abortFlash("Error: Unable to select hyperslab for coordinates memspace\n");
    }*/


    /* read the data */
    if(!*doread || *local_blocks <= 0) {
      ierr = H5Sselect_none(dataspace);
      if(ierr < 0) Driver_abortC("[" FILE_AT_LINE_C "]: hdf5 error.");
      ierr = H5Sselect_none(memspace);
      if(ierr < 0) Driver_abortC("[" FILE_AT_LINE_C "]: hdf5 error.");
    }
#ifdef FLASH_IO_ASYNC_HDF5
    status = H5Dread_async(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                 dxfer_template, unknowns, io_es_id);
#else
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                 dxfer_template, unknowns);
#endif    
               
    if (status < 0){
      printf("Error: Unable to read unknowns\n");
      Driver_abortC("Error: Unable to read unknowns\n");
    }

    H5Pclose(dxfer_template);
    H5Sclose(memspace); 
    H5Sclose(dataspace);
#ifdef FLASH_IO_ASYNC_HDF5
  H5Dclose_async(dataset, io_es_id);
#else
  H5Dclose(dataset);
#endif  
  }
}
