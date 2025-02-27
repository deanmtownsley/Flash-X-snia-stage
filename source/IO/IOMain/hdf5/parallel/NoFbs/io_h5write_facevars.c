#include "mangle_names.h"
#include <hdf5.h>
#include "hdf5_flash.h"
/*#include <mpi.h>*/
#include "Simulation.h"
#include "constants.h"
#include "io_flash.h"
#include "io_h5_attribute.h"

int Driver_abortC(char* message);

/*
   This is based on the methods of writing out unknown arrays.  If the offset
   of a particular axis is not zero, then the read starts one higher than the
   offset passed in for facevar data, to prevent data from separate procs from
   colliding.
   
   This function writes out a single unknown (passed from the checkpoint 
   routine), giving the record a label from the varnames or species
   database 
   
   The dimensions of the unknowns array (nvar, nxb, nyb, nzb, maxblocks)
   are passed through as arguments.  The dataspace (what gets written
   to disk) and the memory space (how the unknowns array is accessed in
   local memory) are defined based on these passed values.  This allows
   use to pass an unk array that has all the guardcells + interior cells
   (as in the checkpointing), or a dummy unk array that has just the 
   interior cells (in this case, nguard would be passed as 0).
*/


void FTOC(io_h5write_facevars)(int* myPE,
                               hid_t* file_identifier,
                               int* globalNXB,       /* num of zones in x dir across whole domain */
                               int* globalNYB,       /* num of zones in y dir across whole domain */
                               int* globalNZB,       /* num of zones in z dir across whole domain */
                               int* nxbOffset,       /* num of zones in x dir to right of 0 where we'll start writing */
                               int* nybOffset,       /* num of zones in y dir to right of 0 where we'll start writing */
                               int* nzbOffset,       /* num of zones in z dir to right of 0 where we'll start writing */
                               int* nxb,             /* num of zones in x dir to be written on this call */
                               int* nyb,             /* num of zones in y dir to be written on this call */
                               int* nzb,             /* num of zones in z dir to be written on this call */
                               double* unknowns,     /* [1][NZB][NYB][NXB] */
			       double* varMin,       /* minimum value of unknown */
			       double* varMax,       /* maximum value of unknown */
                               char record_label[5]) /* add char-null termination */
{

  hid_t dataspace, dataset, memspace, dxfer_template, dataset_plist;

  herr_t status, ierr;

  int rank;
  hsize_t dimens_4d[4], dimens_5d[5];

  hsize_t start_4d[4];
  hsize_t stride_4d[4], count_4d[4];

#ifdef CHUNK
  hsize_t dimens_chunk[4];
#endif

  char record_label_new[5];

  int total_blocks = 1; /* not necessary for nofbs, but keeps consistency for fidlr routines */

  const char minLabel[] = "minimum";
  const char maxLabel[] = "maximum";
  const int attType = IO_FLASH_DOUBLE;
  const int dims = 1;
  const int diskSize[] = {1};

  /*printf("varMin = %e\nvarMax = %e\n", *varMin, *varMax);*/
  /*  printf("myPE: %d\nglobalNXB: %d\nglobalNYB: %d\nglobalNZB: %d\nnxbOffset: %d\nnybOffset: %d\nnzbOffset%d\nnxb: %d\nnyb: %d\nnzb: %d\n",
   *myPE, *globalNXB, *globalNYB, *globalNZB, *nxbOffset, *nybOffset, *nzbOffset, *nxb, *nyb, *nzb);*/



  /* if we have an offset, skip writing the first element.  That has already 
     been written!*/
  /*but this is the wrong way to do it*/
  /*if(nxbOffset > 0)
    *nxbOffset++;
   if(nybOffset > 0)
    *nybOffset++;
  if(nzbOffset > 0)
    *nzbOffset++;
  */
  /* 
     the variable names are 4 characters long -- copy this into 
     record_label_new, the 5th character is for the \0 termination 
  */
  strncpy(record_label_new, record_label,4);
  *(record_label_new + 4) = '\0';

  /* set the dimensions of the dataset */
  rank = 4;
  dimens_4d[0] = total_blocks;
  dimens_4d[1] = *globalNZB;
  dimens_4d[2] = *globalNYB;
  dimens_4d[3] = *globalNXB;

  dataspace = H5Screate_simple(rank, dimens_4d, NULL);
  if(dataspace < 0){
    printf("io_h5write_unknowns: dataspace error");
    Driver_abortC("io_h5write_unknowns: dataspace error");
  }

  dataset_plist = H5Pcreate(H5P_DATASET_CREATE);
  if(dataset_plist < 0){
    printf("io_h5write_unknowns: dataset_plist error");
    Driver_abortC("io_h5write_unknowns: dataset_plist error");
  }

    
  /* create a parallel hdf5 dataset */
  dataset = H5Dcreate(*file_identifier, record_label_new,
	              H5T_NATIVE_DOUBLE, dataspace, dataset_plist); 
  if(dataset < 0) {
    Driver_abortC("dataset Error: H5Dcreate io_h5write_unk\n");
  }    


  /*change the offset here to avoid writing on the same memory locations at the same time, io collision reduction*/
  start_4d[0] = 0;
  if(nzb > nyb && nzb > nxb && nzbOffset != 0)
    start_4d[1] = *nzbOffset + 1; /*number of z zones to the 'left' of myPE -- do this by copying functionality in Particles_getOffset */
  else
    start_4d[1] = *nzbOffset;
  if(nyb > nzb && nyb > nxb && nybOffset != 0)
    start_4d[2] = *nybOffset + 1;
  else
    start_4d[2] = *nybOffset;
  if(nxb > nyb && nxb > nzb && nxbOffset != 0)
    start_4d[3] = *nxbOffset + 1;
  else
    start_4d[3] = *nxbOffset;


  stride_4d[0] = total_blocks;
  stride_4d[1] = 1;
  stride_4d[2] = 1;  
  stride_4d[3] = 1;

  /*now really write the local nxb, nyb and nzb */
  /*remember still 1 block per proc, it's just that the blocks can have various, nxb, nyb, nzb */
  /*ensure that we don't go over our bounds.*/
  count_4d[0] = total_blocks;
  if(nzb > nyb && nzb > nxb && nzbOffset != 0)
    count_4d[1] = *nzb - 1;
  else
    count_4d[1] = *nzb;
  if(nyb > nzb && nyb > nxb && nybOffset != 0)
    count_4d[2] = *nyb - 1;
  else
       count_4d[2] = *nyb;
  if(nxb > nyb && nxb > nzb && nxbOffset != 0)
    count_4d[3] = *nxb - 1; 
  else   
    count_4d[3] = *nxb;

  ierr = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start_4d, 
			     stride_4d, count_4d, NULL);

  if(ierr < 0){
     printf("%s\n", "Error: unable to select hyperslab for unknowns dataspace");
     Driver_abortC("Error: unable to select hyperslab for unknowns dataspace");
  }


  rank = 5;
  dimens_5d[0] = total_blocks;
  dimens_5d[1] = *nzb;
  dimens_5d[2] = *nyb;
  dimens_5d[3] = *nxb;
  dimens_5d[4] = 1;


  memspace = H5Screate_simple(rank, dimens_5d, NULL);
  if(memspace < 0){
    printf("io_h5write_unknowns: memspace error");
    Driver_abortC("io_h5write_unknowns: memspace error");
  }

  /* obtain a copy of the file transfer property list */ 
  dxfer_template = H5Pcreate(H5P_DATASET_XFER);


  /* default for now */
  /* default is 'independent' mode */


  /*ierr = H5Pset_dxpl_mpio(dxfer_template, H5FD_MPIO_INDEPENDENT);*/
  /*ierr = H5Pset_preserve(dxfer_template, 0u);*/

 


#ifdef CHUNK
  /* set the layout to chunked */

  ierr = H5Pset_layout(dataset_plist, H5D_CHUNKED);
  
  /* create a chunk a containing 10 blocks worth of data */
  dimens_chunk[0] = 10;
  dimens_chunk[1] = *nzb;
  dimens_chunk[2] = *nyb;
  dimens_chunk[3] = *nxb;
   
  ierr = H5Pset_chunk(dataset_plist, 4, dimens_chunk);
#endif


  /*!!DEV -kda I've removed the max and min parts for now */

  /* write the data */
  status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, 
                dxfer_template, unknowns);

  if(status < 0){
    printf("io_h5write_unknowns: H5Dwrite error");
    Driver_abortC("io_h5write_unknowns: H5Dwrite error");
  }

#ifdef DEBUG_IO
  printf("UNKNOWNS: wrote unknowns, status = %d\n", (int) status);
#endif

  
  H5Pclose(dxfer_template);

  H5Pclose(dataset_plist);

  H5Sclose(memspace); 
  H5Sclose(dataspace);
  H5Dclose(dataset);


  /* Write attributes for unknowns: varMin and varMax.
     Unlike fixed block size mode, we only have to worry about this in
     parallel for now.  Note: The io_h5_attribute_XXX functions open
     the dataset "record_label_new" before creating / writing
     attributes, so "record_label_new" should be closed before the
     call. */
  io_h5_attribute_create(*myPE, (int)*file_identifier, attType, dims,
			 diskSize, record_label_new, minLabel);
  io_h5_attribute_write(*myPE, (int)*file_identifier, attType,
			record_label_new, minLabel, varMin);

  io_h5_attribute_create(*myPE, (int)*file_identifier, attType, dims,
			 diskSize, record_label_new, maxLabel);
  io_h5_attribute_write(*myPE, (int)*file_identifier, attType,
			record_label_new, maxLabel, varMax);
}
