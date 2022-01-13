#include "io_h5_attribute.h"
#include "Simulation.h"

/* WARNING: The caller must null-terminate all strings */
#ifdef FLASH_IO_ASYNC_HDF5
  hid_t io_es_id;
#endif

void io_h5_attribute_create(const int myPE,
			    const hid_t fileID,
			    const int diskType,
			    const int dims,
			    const int diskSize[],
			    const char datasetName[],
			    const char attDatasetName[])
{
  const hid_t hFileID = fileID;
  hsize_t hDiskSize[IO_MAX_DIMS];
  hid_t hDiskType, dsetID, attID, attDataspace;
  herr_t err;
  int trueDims, parallelIO, i;

#if defined(IO_HDF5_PARALLEL)
  parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
  parallelIO = 0;
#else
  parallelIO = 0;
#endif

  if (parallelIO || (!parallelIO && myPE == MASTER_PE)) {
    assert(dims > 0 && dims <= IO_MAX_DIMS);
    for (i=0; i<dims; ++i) {
      hDiskSize[i] = (hsize_t) diskSize[i];
    }

    trueDims = dims;
    if (diskType == IO_FLASH_STRING) {
      hDiskType = io_h5_type_create_string(diskSize[dims-1]);
      if (dims == 1) {
	/* There is only 1 string */
	hDiskSize[0] = 1;
      } else if (dims > 1) {
	trueDims = trueDims - 1;
      }
    } else {
      /* We are using simple primitive types */
      hDiskType = io_h5_type_hid_primitive(diskType);
    }

    /* Open appropriate dataset and then add the attribute */
#ifdef FLASH_IO_ASYNC_HDF5
    dsetID = H5Dopen_async(hFileID, datasetName,H5P_DEFAULT, io_es_id);
#else
    dsetID = H5Dopen(hFileID, datasetName,H5P_DEFAULT);
#endif
    assert(dsetID >= 0);

    attDataspace = H5Screate_simple(trueDims, hDiskSize, NULL);
    assert(attDataspace >= 0);

#ifdef FLASH_IO_ASYNC_HDF5
    attID = H5Acreate_async(dsetID, attDatasetName, hDiskType,
		      attDataspace, H5P_DEFAULT,H5P_DEFAULT, io_es_id);
#else
    attID = H5Acreate(dsetID, attDatasetName, hDiskType,
		      attDataspace, H5P_DEFAULT,H5P_DEFAULT);
#endif

    assert(attID >= 0);


    /* Free allocated space. */
    err = H5Sclose(attDataspace);
    assert(err >= 0);

#ifdef FLASH_IO_ASYNC_HDF5
    err = H5Aclose_async(attID, io_es_id);
#else
    err = H5Aclose(attID);
#endif
    assert(err >= 0);

#ifdef FLASH_IO_ASYNC_HDF5
    err = H5Dclose_async(dsetID, io_es_id);
#else
    err = H5Dclose(dsetID);
#endif
    assert(err >= 0);

    if (diskType == IO_FLASH_STRING) {
      io_h5_type_free_string(hDiskType);
    }
  }
}


/* We pass the memory datatype so that the HDF5 library can
   convert between double in memory and float in file */
void io_h5_attribute_write(const int myPE,
			   const hid_t fileID,
			   const int memType,
			   const char datasetName[],
			   const char attDatasetName[],
			   const void * const pData)
{
  const hid_t hFileID = fileID;
  hid_t hMemType, dsetID, attID;
  herr_t err;
  int parallelIO;

#if defined(IO_HDF5_PARALLEL)
  parallelIO = 1;
#elif defined(IO_HDF5_SERIAL)
  parallelIO = 0;
#else
  parallelIO = 0;
#endif

  if (parallelIO || (!parallelIO && myPE == MASTER_PE)) {

    /* Open appropriate dataset and then attribute */
#ifdef FLASH_IO_ASYNC_HDF5
    dsetID = H5Dopen_async(hFileID, datasetName,H5P_DEFAULT, io_es_id);
#else
    dsetID = H5Dopen(hFileID, datasetName,H5P_DEFAULT);
#endif

    assert(dsetID >= 0);


#if (H5_VERS_MAJOR == 1 && H5_VERS_MINOR < 8)
    /* This function is now deprecated by HDF5 */
    attID = H5Aopen_name(dsetID, attDatasetName);
#else
#ifdef FLASH_IO_ASYNC_HDF5
    attID = H5Aopen_async(dsetID, attDatasetName, H5P_DEFAULT, io_es_id);
#else
    attID = H5Aopen(dsetID, attDatasetName, H5P_DEFAULT);
#endif
#endif
    assert(attID >= 0);


    if (memType == IO_FLASH_STRING) {
      hMemType = H5Aget_type(attID);
      assert(hMemType >= 0);
    } else {
      hMemType = io_h5_type_hid_primitive(memType);
    }

    /* Write the attribute */
#ifdef FLASH_IO_ASYNC_HDF5
    err = H5Awrite_async(attID, hMemType, pData, io_es_id);
#else
    err = H5Awrite(attID, hMemType, pData);
#endif

    assert(err >= 0);


    if (memType == IO_FLASH_STRING) {
      io_h5_type_free_string(hMemType);
    }

#ifdef FLASH_IO_ASYNC_HDF5
    err = H5Aclose_async(attID, io_es_id);
#else
    err = H5Aclose(attID);
#endif

    assert(err >= 0);
#ifdef FLASH_IO_ASYNC_HDF5
    err = H5Dclose_async(dsetID, io_es_id);
#else
    err = H5Dclose(dsetID);
#endif
    assert(err >= 0);
  }
}
