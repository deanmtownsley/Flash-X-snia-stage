#include <unistd.h>
#include <stdlib.h>
#include "options.h"

int parse_cmdline(options_t *opts, int argc, char **argv)
{
  extern char *optarg;
  extern int optind, optopt;
  char *file1, *file2, *pwd;
  int i, j, pwdLen, filenameLen;
  int argError = 0, c;

  /* initialize opts */
  opts->ignore = 0;
  opts->pComp = 0;
  opts->verbose = 0;
  opts->reorder = 0;
  opts->benchmark = 'B';
  opts->norm_order = 1;
  opts->mesh_tol = 0.;
  opts->err_tol = 0.;
  opts->perr_tol = 0.;
  opts->pmatch_tol = 0.;
  opts->extra_varnames = 0;
  opts->n_extra_varnames = 0;
  opts->ignorableVarnames = 0;
  opts->n_ignorableVarnames = 0;
  opts->pCompVarnames = 0;
  opts->n_pCompVarnames = 0;
  opts->gridVarSelfDiscovery = 0;

  /* Argument parsing */
  while ((c = getopt(argc, argv, "ahvirb:L:t:s:I:D:e:p:P:")) != -1)
  {
    switch (c)
    {
    case 'a':
      opts->gridVarSelfDiscovery = 1;
      break;
    case 'h':
      argError = 1;
      break;
    case 'v':
      opts->verbose = 1;
      break;
    case 'i':
      opts->ignore = 1;
      break;
    case 'r':
      opts->reorder = 1;
      break;
    case 'L':
      opts->norm_order = atoi(optarg);
      break;
    case 't':
      opts->mesh_tol = atof(optarg);
      break;
    case 'e':
      opts->err_tol = atof(optarg);
      break;
    case 'p':
      opts->perr_tol = atof(optarg);
      break;
    case 'P':
      opts->pmatch_tol = atof(optarg);
      break;
    case 'b':
      opts->benchmark = optarg[0];
      if (opts->benchmark == 'a')
      {
        opts->benchmark = 'A';
      }
      if (opts->benchmark == 'b')
      {
        opts->benchmark = 'B';
      }
      break;
    case 's':
      if (strlen(optarg) > 4)
        return 1;
      if (!opts->extra_varnames)
        opts->extra_varnames = (char **)malloc(sizeof(char *));
      else
        opts->extra_varnames = (char **)realloc(opts->extra_varnames, sizeof(char *) * (opts->n_extra_varnames + 1));
      opts->extra_varnames[opts->n_extra_varnames] = (char *)malloc(5);
      strncpy(opts->extra_varnames[opts->n_extra_varnames], optarg, 5);
      opts->extra_varnames[opts->n_extra_varnames][4] = (char)0;
      for (i = 1; i < 4; i++)
      {
        if (opts->extra_varnames[opts->n_extra_varnames][i] == (char)0)
          opts->extra_varnames[opts->n_extra_varnames][i] = ' ';
      }
      opts->n_extra_varnames++;
      break;
    case 'I':
      if (strlen(optarg) > 4)
        return 1;
      if (!opts->ignorableVarnames)
        opts->ignorableVarnames = (char **)malloc(sizeof(char *));
      else
        opts->ignorableVarnames = (char **)realloc(opts->ignorableVarnames, sizeof(char *) * (opts->n_ignorableVarnames + 1));
      opts->ignorableVarnames[opts->n_ignorableVarnames] = (char *)malloc(5);
      strncpy(opts->ignorableVarnames[opts->n_ignorableVarnames], optarg, 5);
      opts->ignorableVarnames[opts->n_ignorableVarnames][4] = (char)0;
      for (i = 1; i < 4; i++)
      {
        if (opts->ignorableVarnames[opts->n_ignorableVarnames][i] == (char)0)
          opts->ignorableVarnames[opts->n_ignorableVarnames][i] = ' ';
      }
      opts->n_ignorableVarnames++;
      break;
    case 'D':
      opts->pComp = 1;
      if (strlen(optarg) > 4)
        return 1;
      if (!opts->pCompVarnames)
        opts->pCompVarnames = (char **)malloc(sizeof(char *));
      else
        opts->pCompVarnames = (char **)realloc(opts->pCompVarnames, sizeof(char *) * (opts->n_pCompVarnames + 1));
      opts->pCompVarnames[opts->n_pCompVarnames] = (char *)malloc(5);
      strncpy(opts->pCompVarnames[opts->n_pCompVarnames], optarg, 5);
      opts->pCompVarnames[opts->n_pCompVarnames][4] = (char)0;
      for (i = 1; i < 4; i++)
      {
        if (opts->pCompVarnames[opts->n_pCompVarnames][i] == (char)0)
          opts->pCompVarnames[opts->n_pCompVarnames][i] = ' ';
      }
      opts->n_pCompVarnames++;
      break;
    case '?':
      /* do nothing - this is a change to prevent sfocu from
         crashing when it tries to interpret options meant for
         mpirun or poe. I have left the original 'return 1' in
         this commented-out block
      return 1;
      */
      break;
    }
  }

  if ((argc - optind) != 2)
    argError = 1;
  else
  {
    pwd = getenv("PWD");
    pwdLen = strlen(pwd);
    if (argv[optind][0] == '/')
    {
      file1 = argv[optind];
    }
    else
    {
      /* we have a relative path, so prepend the pwd */
      filenameLen = strlen(argv[optind]);
      file1 = (char *)malloc(sizeof(char) * (pwdLen + filenameLen + 2));
      for (i = 0; i < pwdLen; i++)
      {
        file1[i] = pwd[i];
      }
      file1[i] = '/';
      i++;
      for (j = 0; j < filenameLen; j++)
      {
        file1[i] = argv[optind][j];
        i++;
      }
      file1[i] = 0;
    }
    if (argv[optind + 1][0] == '/')
    {
      file2 = argv[optind + 1];
    }
    else
    {
      /* we have a relative path, so prepend the pwd */
      filenameLen = strlen(argv[optind + 1]);
      file2 = (char *)malloc(sizeof(char) * (pwdLen + filenameLen + 2));
      for (i = 0; i < pwdLen; i++)
      {
        file2[i] = pwd[i];
      }
      file2[i] = '/';
      i++;
      for (j = 0; j < filenameLen; j++)
      {
        file2[i] = argv[optind + 1][j];
        i++;
      }
      file2[i] = 0;
    }
    opts->file1 = file1;
    opts->file2 = file2;
  }

  return argError;
}
