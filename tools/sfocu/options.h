#ifndef _OPTIONS_H
#define _OPTIONS_H

#include <string.h>
typedef struct options_t
{
  int verbose;
  int ignore;
  int pComp;
  int norm_order;
  int reorder;
  double mesh_tol;
  double err_tol;
  double perr_tol;
  double pmatch_tol;
  char *file1, *file2;
  char **extra_varnames;
  int n_extra_varnames;
  char **ignorableVarnames;
  int n_ignorableVarnames;
  char **pCompVarnames;
  int n_pCompVarnames;
  int gridVarSelfDiscovery;
  char benchmark;
} options_t;

int parse_cmdline(options_t *opts, int argc, char **argv);

#endif /* _OPTIONS_H */
