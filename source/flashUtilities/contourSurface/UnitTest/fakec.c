
#include <stdio.h>
#include <stdlib.h>

int Driver_abortC(char* message){
  fprintf(stderr, message);
  exit(1);

  return 0;
}

