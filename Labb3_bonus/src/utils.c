#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "perc.h"


int write_data(Sys_Params *params, int n_percs, char *fname) {
    int N = params->N;

    char path[FNAMESIZE + 5] = "data/";
    strcat(path, fname);
    
    FILE *fp;

    if (fp = fopen(path, "a")) {
        fprintf(fp, "%d %d\n", n_percs, N);
        fclose(fp);
    }

    return 1;
}

