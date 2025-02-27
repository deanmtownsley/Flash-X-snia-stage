#ifndef _FLASH_READER_H
#define _FLASH_READER_H

#define FR_MDIM 3
#define FR_VAR_STRING_SIZE 4
#define FR_LEAF_NODE 1
#define FR_MAXVARS 500 
#define FR_MAXPARTPROPS 100
#include "options.h"

/* Supported formats */
#define FR_HDF5   1
#define FR_HDF5_PMESH 4
#define FR_HDF5_CHOMBO 5
#define FR_HDF4 2
#define FR_NCDF 3

#define FR_PART_PROP_STRING_SIZE 24
#define MAX_STRING_LENGTH 80

/*Define our variable types*/
#define UNK 1
#define SCRATCH 2
#define FACEX 3
#define FACEY 4
#define FACEZ 5

#define XDIM 3
#define YDIM 2
#define ZDIM 1
#define NUMBLOCKS 0

/* this represents a collection of particles and all
   properties */
typedef struct FR_ParticlesAllProps {
  int numParticles;      /* number of particles in this collection */
  int numIntProps;       /* number of int props for each particle */
  int numRealProps;
  int *intProps;         /* big contiguous array of int prop values */
  double *realProps;
  char **intPropsNames;  /* names the props in the contiguous values array */
  char **realPropsNames;
} FR_ParticlesAllProps;

typedef struct FR_File{
  int format;
  char filename[1024];
  int nblocks;
  /* ncells are only valid for files where all blocks are same size (not chombo) */
  int ncells_vec[FR_MDIM];
  int ncells; /* # cells in a block: product of ncells_vec's components*/
  int dim;
  int *nodetype;
  double *coord;
  double *size;
  double *bbox;
  int *lref;
  int nvar;
  char varnames[FR_MAXVARS][FR_VAR_STRING_SIZE+1];
  int vartypes[FR_MAXVARS];
  int totalparticles;
  int numRealPartProps;
  char realPartPropNames[FR_MAXVARS][FR_PART_PROP_STRING_SIZE+1];
  int numIntPartProps;
  char intPartPropNames[FR_MAXVARS][FR_PART_PROP_STRING_SIZE+1];

  int *blocksPerRefineLevel; /*For chombo hdf5 format */	
  double *levelSize;
  int nLevels;
  

  void *handle; /* This points to a handle appropriate for the file's format */
} FR_File;

typedef struct FR_Block{
  double *data; 
  int size[FR_MDIM];
  /* 
  int dim1, dim2, dim3;
     Get data in cell (nxb, nyb, nzb) through data[nxb*dim1+nyb*dim2+nzb*dim3]
     Inefficient but simple. Do this later since for now we can do the 
     comparisons without knowing how the data array is arranged.
  */
} FR_Block;

extern FR_File *FR_open(char *filename, options_t opts);
extern void FR_close(FR_File *);
extern FR_Block *FR_GetBlock(FR_File *file, char *var, int block_no);
extern FR_Block *FR_GetBlock_HDF5_Chombo(FR_File *file, int varPos, int block_no);
extern FR_ParticlesAllProps *FR_GetParticlesAllProps(FR_File *file, int startingAt, int atMost, int *numberGotten);
extern void FR_DeleteBlock(FR_Block *);
extern void FR_DeleteParticlesAllProps(FR_ParticlesAllProps *);

#endif /* _FLASH_READER_H */
