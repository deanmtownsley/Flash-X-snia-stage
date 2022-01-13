#include "options.h"
#include "flash_reader.h"

#define NO_BLOCK -1 /* Not one of paramesh's nodetypes, I think */

typedef struct particleReport
{
  int particlesCompared; /* number of particles we compare */
  int extraParticlesA;   /* particles in A but not B */
  int extraParticlesB;
  int numRealPartProps;         /* number of properties we compare */
  int fail;                     /* if false, the files are "identical" */
  int badParticles[FR_MAXPARTPROPS]; /* number of particles where a prop differs by > perr_tol */
  double sumA[FR_MAXPARTPROPS], sumB[FR_MAXPARTPROPS];
  double maxA[FR_MAXPARTPROPS], maxB[FR_MAXPARTPROPS];
  double minA[FR_MAXPARTPROPS], minB[FR_MAXPARTPROPS];
  double maxError[FR_MAXPARTPROPS]; /* max of error norm */
  double maxErrorVal[2][FR_MAXPARTPROPS];
  double minError[FR_MAXPARTPROPS];
  double minErrorVal[2][FR_MAXPARTPROPS]; /*we calculate but don't report this*/
  double absError[FR_MAXPARTPROPS];       /* sup |a-b| */
  double absErrorVal[2][FR_MAXPARTPROPS];
  double magError[FR_MAXPARTPROPS];                                      /* sup |a-b| / max( sup |a|, sup |b|, 1e-99)*/
  char realPartPropNames[FR_MAXPARTPROPS][FR_PART_PROP_STRING_SIZE + 1]; /* properties we compare */
  int ipropA[FR_MAXPARTPROPS], ipropB[FR_MAXPARTPROPS];                       /* maps from report props to actual file props */
} particleReport;

void compareparticles(options_t *opts, FR_File *A, FR_File *B, FR_ParticlesAllProps *pAllA, FR_ParticlesAllProps *pAllB, int *pAtoB, struct particleReport *particleReport, int particlesGotten);
