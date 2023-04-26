#include "compareparticles.h"
#include "flash_reader.h"
#include "options.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

/* find file A and file B's <prop> index */
int find_particle_prop_index(FR_ParticlesAllProps *pAll, char *propName) {
  int retIndex = -1, propIndex;
  for (propIndex = 0; propIndex < pAll->numRealProps; propIndex++) {
    if (strncmp(pAll->realPropsNames[propIndex], propName, strlen(propName)) == 0) {
      retIndex = propIndex;
      break;
    }
  }
  return retIndex;
}

/* norm(a, b)

   Same norm as FOCU; if you change this change the output message below as
   well. This is called only when a!=b 
*/
double norm(double a, double b) {
  return fabs(2 * (a - b)) / (fabs(a + b) > 1e-99 ? fabs(a + b) : 1e-99);
}

void compareparticles(options_t *opts, FR_File *A, FR_File *B, FR_ParticlesAllProps *pAllA, FR_ParticlesAllProps *pAllB, int *pAtoB, struct particleReport *particleReport, int particlesGotten) {
  int i, j, iprop, partIndexA, partIndexB;
  int failed = 0;

  int tagIndexA = -1, tagIndexB = -1, cpuIndexA = -1, cpuIndexB = -1;
  int blkIndexA = -1, procIndexA = -1, tempIndexA = -1, tempIndexB = -1;
  int posxIndexA = -1, posxIndexB = -1, posyIndexA = -1, posyIndexB = -1, poszIndexA = -1, poszIndexB = -1;
  int velxIndexA = -1, velxIndexB = -1, velyIndexA = -1, velyIndexB = -1, velzIndexA = -1, velzIndexB = -1;
  int propIndexA[FR_MAXVARS], propIndexB[FR_MAXVARS];
  double pPropValueA[FR_MAXVARS], pPropValueB[FR_MAXVARS];

  bool found = false;
  int nPVars = 0;
  double error = 0.0;
  double onePropDiff;

  double propValueA, propValueB, posxValueA, posxValueB, posyValueA, posyValueB, poszValueA, poszValueB;

  /* check for particle cpu prop */
  tagIndexA = find_particle_prop_index(pAllA, "tag");
  tagIndexB = find_particle_prop_index(pAllB, "tag");

  cpuIndexA = find_particle_prop_index(pAllA, "cpu");
  cpuIndexB = find_particle_prop_index(pAllB, "cpu");

  /*find the particle blk index 
        we don't want to compare this property, because simulations on different
        number of procs will have particles with different blk properties even
        though the particles are in the same physical domain. */

  blkIndexA = find_particle_prop_index(pAllA, "blk");
  procIndexA = find_particle_prop_index(pAllA, "proc");

  /* opts->pComp is zero when -D options isn't specified
   cpuIndex > -1 when particles have a cpu tag. */
  if ((cpuIndexA > -1 || cpuIndexB > -1) && opts->pComp == 0) {
    if (opts->verbose)
      printf("Searching particles by comparing pos/vel {x,y,z}\n");
    posxIndexA = find_particle_prop_index(pAllA, "posx");
    posxIndexB = find_particle_prop_index(pAllB, "posx");

    posyIndexA = find_particle_prop_index(pAllA, "posy");
    posyIndexB = find_particle_prop_index(pAllB, "posy");

    poszIndexA = find_particle_prop_index(pAllA, "posz");
    poszIndexB = find_particle_prop_index(pAllB, "posz");

    velxIndexA = find_particle_prop_index(pAllA, "velx");
    velxIndexB = find_particle_prop_index(pAllB, "velx");

    velyIndexA = find_particle_prop_index(pAllA, "vely");
    velyIndexB = find_particle_prop_index(pAllB, "vely");

    velzIndexA = find_particle_prop_index(pAllA, "velz");
    velzIndexB = find_particle_prop_index(pAllB, "velz");
    /* set the counter for loop comparison later 3 vel 3 pos*/
    nPVars = 6;
    /* populate index for both files*/
    propIndexA[0] = posxIndexA;
    propIndexB[0] = posxIndexB;
    propIndexA[1] = posyIndexA;
    propIndexB[1] = posyIndexB;
    propIndexA[2] = poszIndexA;
    propIndexB[2] = poszIndexB;
    propIndexA[3] = velxIndexA;
    propIndexB[3] = velxIndexB;
    propIndexA[4] = velyIndexA;
    propIndexB[4] = velyIndexB;
    propIndexA[5] = velzIndexA;
    propIndexB[5] = velzIndexB;
  }

  /* Let's get prop Index based on -D option specified */
  for (i = 0; i < opts->n_pCompVarnames; i++) {
    nPVars = opts->n_pCompVarnames;
    tempIndexA = find_particle_prop_index(pAllA, opts->pCompVarnames[i]);
    tempIndexB = find_particle_prop_index(pAllB, opts->pCompVarnames[i]);
    propIndexA[i] = tempIndexA;
    propIndexB[i] = tempIndexB;
  }

  /* start of particle table-building loop */
      /* compute sorted permutation of B particle tags for fast lookup */
      /* used for if 0 - particle lookup by sorted tag case. */
#if 0      
      int *pBsort = sort_particle_tags(
        pAllB->realProps, B->totalparticles, B->numRealPartProps, tagIndexB
      );
#endif

  for (partIndexA = 0; partIndexA < A->totalparticles; partIndexA++) {
#if 0 /* use particle lookup by sorted tag */
        partIndexB = find_particle_tag(
          pAllA->realProps[partIndexA*pAllA->numRealProps + tagIndexA],
          pAllB->realProps, pBsort, B->totalparticles, B->numRealPartProps, tagIndexB
        );
        if(partIndexB >= 0 && 
          pAllA->realProps[partIndexA*pAllA->numRealProps + tagIndexA]
          != pAllB->realProps[partIndexB*pAllB->numRealProps + tagIndexB]) {
          printf("BUG in particle sort!\n");
        }
        pAtoB[partIndexA] = partIndexB;
        failed = partIndexB < 0;
#else /* do not use particle tag sorting */
    /* this is common */
    if (cpuIndexA == -1 && cpuIndexB == -1) {
      if (pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA] ==
          pAllB->realProps[(partIndexA * pAllB->numRealProps) + tagIndexB]) {
        pAtoB[partIndexA] = partIndexA;
        continue;
      }
      failed = 1;
      for (partIndexB = 0; partIndexB < B->totalparticles; partIndexB++) {
        if (pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA] ==
            pAllB->realProps[(partIndexB * pAllB->numRealProps) + tagIndexB]) {
          pAtoB[partIndexA] = partIndexB;
          failed = 0;
          break;
        }
      }
    } else if (cpuIndexA > -1 || cpuIndexB > -1) { /* we have cpu property, particles equality is different here */
      found = false;
      for (j = 0; j < nPVars; j++) {
	pPropValueA[j] = pAllA->realProps[(partIndexA * pAllA->numRealProps) + propIndexA[j]];
	pPropValueB[j] = pAllB->realProps[(partIndexA * pAllB->numRealProps) + propIndexB[j]];

	if (pPropValueA[j] == pPropValueB[j]) { /* exact match */
	  found = true;
	} else if (opts->pmatch_tol > 0.0) { /* try tolerance */
	  onePropDiff = fabs(pPropValueA[j] - pPropValueB[j]);
	  if (onePropDiff <= opts->pmatch_tol)
	    found = true;
	  else {
	    found = false;
	    break;
	  }
	} else {
	  found = false;
	  break;
	}
      }
      if (found) {
        pAtoB[partIndexA] = partIndexA;
	continue;
      }
      failed = 1;
      for (partIndexB = 0; partIndexB < B->totalparticles; partIndexB++) {
        for (j = 0; j < nPVars; j++) {
          pPropValueA[j] = pAllA->realProps[(partIndexA * pAllA->numRealProps) + propIndexA[j]];
          pPropValueB[j] = pAllB->realProps[(partIndexB * pAllB->numRealProps) + propIndexB[j]];

          if (pPropValueA[j] == pPropValueB[j]) { /* exact match */
            found = true;
	  } else if (opts->pmatch_tol > 0.0) { /* try tolerance */
	    onePropDiff = fabs(pPropValueA[j] - pPropValueB[j]);
	    if (onePropDiff <= opts->pmatch_tol)
	      found = true;
	    else {
	      found = false;
	      break;
	    }
          } else {
            found = false;
            break;
          }
        }
        failed = 1;
        if (found == true) {
          if (opts->verbose) {
            printf("A: velx %f, posx, %f\n", pPropValueA[0], pPropValueA[3]);
            printf("B: velx %f, posx, %f\n", pPropValueB[0], pPropValueB[3]);
            printf("A%d, -> B%d\n", partIndexA, partIndexB);
          }
          pAtoB[partIndexA] = partIndexB;
          failed = 0;
          break;
        }
      }
    } /* end else cpuIndex block */
#endif
    if (failed) {
      pAtoB[partIndexA] = -1;
      particleReport->extraParticlesA++;
      particleReport->fail = 1;

      /*      if(opts->verbose) */
      printf("No match for particle %d (tag=%16.16g) in %s\n",
             partIndexA,
             pAllA->realProps[(partIndexA * pAllA->numRealProps) + tagIndexA],
             A->filename);
      /*exit(1);*/
    }
  } /* end of particle table-building loop */

#if 0
    free(pBsort);        
#endif


  /* Now Check that particle real properties are the same */
  particleReport->particlesCompared = 0;
  for (i = 0; i < particlesGotten; i++) {
    if (pAtoB[i] != NO_BLOCK) {
      for (iprop = 0; iprop < particleReport->numRealPartProps; iprop++) {
	int ipropA = particleReport->ipropA[iprop];
        /*AMReX case or case with cpu id */
        if (cpuIndexA > -1 || cpuIndexA > -1) {
          /* don't check for tag or cpu index */
          if (ipropA == tagIndexA || ipropA == cpuIndexA) {
            continue;
          }
        }
        if (ipropA != blkIndexA && ipropA != procIndexA) {
          propValueA = pAllA->realProps[i * A->numRealPartProps + ipropA];
          propValueB = pAllB->realProps[pAtoB[i] * B->numRealPartProps + particleReport->ipropB[iprop]];

          if (particleReport->particlesCompared == 0) {
            particleReport->sumA[iprop] = 0;
            particleReport->maxA[iprop] = propValueA;
            particleReport->minA[iprop] = propValueA;

            particleReport->sumB[iprop] = 0;
            particleReport->maxB[iprop] = propValueB;
            particleReport->minB[iprop] = propValueB;

            if (propValueA != propValueB)
              particleReport->minError[iprop] = norm(propValueA, propValueB);
            else
              particleReport->minError[iprop] = 0;

            particleReport->badParticles[iprop] = 0;
            particleReport->maxError[iprop] = 0.0;
            particleReport->absError[iprop] = 0.0;
          }

          if (propValueA != propValueB) {

            particleReport->badParticles[iprop]++;

            if (opts->verbose)
              printf(" %d %s\n", i, particleReport->realPartPropNames[iprop]);
            error = norm(propValueA, propValueB);

            if (error > particleReport->maxError[iprop]) {
              particleReport->maxError[iprop] = error;
              particleReport->maxErrorVal[0][iprop] = propValueA;
              particleReport->maxErrorVal[1][iprop] = propValueB;
            }
          }

          if (error < particleReport->minError[iprop]) {
            particleReport->minError[iprop] = error;
            particleReport->minErrorVal[0][iprop] = propValueA;
            particleReport->minErrorVal[1][iprop] = propValueB;
          }

          error = fabs(propValueA - propValueB);
          if (error > particleReport->absError[iprop]) {
            particleReport->absError[iprop] = error;
            particleReport->absErrorVal[0][iprop] = propValueA;
            particleReport->absErrorVal[1][iprop] = propValueB;
          }

          particleReport->sumA[iprop] += propValueA;
          particleReport->maxA[iprop] = (propValueA > particleReport->maxA[iprop]) ? propValueA : particleReport->maxA[iprop];
          particleReport->minA[iprop] = (propValueA < particleReport->minA[iprop]) ? propValueA : particleReport->minA[iprop];

          particleReport->sumB[iprop] += propValueB;
          particleReport->maxB[iprop] = (propValueB > particleReport->maxB[iprop]) ? propValueB : particleReport->maxB[iprop];
          particleReport->minB[iprop] = (propValueB < particleReport->minB[iprop]) ? propValueB : particleReport->minB[iprop];

        } /*if iprop != blk prop */
      }

      particleReport->particlesCompared++;
    } /* if pAtoB[i] is valid */
  }
}
