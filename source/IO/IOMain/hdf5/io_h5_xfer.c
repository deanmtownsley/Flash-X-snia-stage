#include "io_h5_xfer.h"
#include "mangle_names.h"

#ifdef FLASH_IO_ASYNC_HDF5
  extern hid_t io_es_id;
#endif

int io_h5_xfer(const int myPE, const hid_t fileID, const int xferType,
	       const hid_t hXferList, const char datasetName[],
	       const hid_t hMemType, const hsize_t hMemSize[],
	       const hsize_t hMemStart[], const hsize_t hMemCount[],
	       const hsize_t hDiskStart[], const hsize_t hDiskCount[],
         const int dims, const int numFileBlks[], void * pData)
{
  herr_t (*old_func)(void*) = NULL;
  void *old_client_data = NULL;
  hsize_t hDummyMemSize[IO_MAX_DIMS];
  hid_t dataspace, memspace, dataset=0;
  herr_t err;
  int i, zeroSize = 0;
  int ierr;


  for (i=0; i<dims; ++i) {
    /* This function only deals with primitive types so following assert
       must be true */
    assert (hMemCount[i] == hDiskCount[i]);
    if (hMemCount[i] == 0 || hDiskCount[i] == 0) zeroSize = 1;
  }

  // search for the dataset
  if (xferType == IO_READ_XFER) {
    ierr = H5Lexists(fileID, datasetName, H5P_DEFAULT);
    
    if (ierr < 0) {
      dataset = -1;
    }
  }

if (dataset == 0){
  /* open the dataset*/
#ifdef FLASH_IO_ASYNC_HDF5
  dataset = H5Dopen_async(fileID, datasetName, H5P_DEFAULT, io_es_id);
#else
  dataset = H5Dopen(fileID, datasetName, H5P_DEFAULT);
#endif  
}

  if (dataset < 0) {

    if (xferType == IO_READ_XFER) {
      ierr = -1;
      if (myPE == MASTER_PE) {
	printf(" [%s]: Skipping missing dataset '%s'.\n",
	       __FILE__, datasetName);
      }
    } else {
      printf("[%s]: Processor %d failed during H5Dopen on dataset %s.\n",
	     __FILE__, myPE, datasetName);
      ierr = Driver_abortC("Error! H5Dopen failed");
    }

  } else {

    // dataspace = H5Dget_space(dataset);
    // assert (dataspace >= 0);


    for (i=0; i<dims; ++i) {
      if (zeroSize) {
	hDummyMemSize[i] = 1; /* Must be > 0 to satisfy H5Screate_simple */
      } else {
	hDummyMemSize[i] = hMemSize[i];
      }
    }

    memspace = H5Screate_simple(dims, hDummyMemSize, NULL);
    assert (memspace >= 0);

    // Convert to hsize_t
    hsize_t hnumFileBlks[IO_MAX_DIMS];
    for (i=0; i<dims; ++i) {
        hnumFileBlks[i] = (hsize_t)(numFileBlks[i]);
    }

    dataspace = H5Screate_simple(dims, hnumFileBlks, NULL);
    assert (dataspace >= 0);

    if (zeroSize) {
      /* Select nothing if myPE is not contributing any data to disk. */
      err = H5Sselect_none(dataspace);
      assert (err >= 0);
      err = H5Sselect_none(memspace);
      assert (err >= 0);
    } else {
      /* Define myPE's portion in global data */
      err = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, hDiskStart,
				NULL, hDiskCount, NULL);
      assert (err >= 0);

      /* Define how myPE's data is stored in memory */
      err = H5Sselect_hyperslab(memspace, H5S_SELECT_SET,
				hMemStart, NULL, hMemCount, NULL);
      assert (err >= 0);
    }


    if (xferType == IO_READ_XFER) {
      /* Read data from disk hyperslab to memory */
#ifdef FLASH_IO_ASYNC_HDF5
      err = H5Dread_async(dataset, hMemType, memspace, dataspace,
        hXferList, pData, io_es_id);
#else
      err = H5Dread(dataset, hMemType, memspace, dataspace,
		    hXferList, pData);
#endif      

    } else if (xferType == IO_WRITE_XFER) {
      /* Write data from memory to disk hyperslab */
#ifdef FLASH_IO_ASYNC_HDF5
      err = H5Dwrite_async(dataset, hMemType, memspace, dataspace,
	            hXferList, pData, io_es_id);
#else
      err = H5Dwrite(dataset,hMemType, memspace, dataspace,
		    hXferList, pData);
#endif         
    }
    assert (err >= 0);


    err = H5Sclose(dataspace);
    assert (err >= 0);
#ifdef FLASH_IO_ASYNC_HDF5
    err = H5Dclose_async(dataset, io_es_id);
#else
    err = H5Dclose(dataset);
#endif
    assert (err >= 0);

    err = H5Sclose(memspace);
    assert (err >= 0);

    ierr = 0;
  }

  return ierr;
}
