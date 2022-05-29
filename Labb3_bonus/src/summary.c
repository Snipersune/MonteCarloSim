#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "perc.h"

#define NVAR 2

int read_data(Sys_Params *params, double *vsum, FILE *stream) {
    int nlines = 0;
    double v[NVAR];

    for (int i = 0; i < NVAR; i++) {
        vsum[i] = 0;
    }

    while (fscanf(stream, "%lf %lf\n", &v[0], &v[1]) != EOF) {
        for (int i = 0; i < NVAR; i++) {
            vsum[i] += v[i];
        }
        nlines++;
    }

    return nlines;
}

void result(Sys_Params *params, double *vsum) {
  
    double N_perc, N_samps;
    N_perc = vsum[0]; 
    N_samps = vsum[1];

    double P_mean, P2_mean, P_err;

    // Get expected P mean and expected P squared mean
    P_mean = N_perc/N_samps;
    P2_mean = P_mean;           // Will be same (if x=1 => x = x^2)


    // Standard deviation in estimated mean of P
    P_err = sqrt(1/(N_samps-1.0)*(P2_mean - P_mean*P_mean));

    fprintf(stdout, "%6e %6e %6e\n", params->p, P_mean, P_err);
}



int main(int argc, char *argv[]) {

    char fname[FNAMESIZE];
    int nblock;

    double vsum[NVAR];
    
    Sys_Params *params = malloc(sizeof(Sys_Params));

    while (fgets(fname, FNAMESIZE, stdin)) {

        // Change the trailing '\n' into a null character.
        char *s = strchr(fname, '\n');
        if (s - fname < FNAMESIZE)
            *s = '\0';

        // It is easier to use files with formatted data (written with fprintf).
        FILE *stream;
        stream = fopen(fname, "r");
        if (!stream) {
            fprintf(stderr, "*** Can't open file %s.\n", fname);
            continue;
        }
        fscanf(stream, "Lx=%d Ly=%d Lz=%d p=%le\n", &params->Lx, &params->Ly, &params->Lz, &params->p);

        if (read_data(params, vsum, stream))
            result(params, vsum);

        fclose(stream);
    }

  free(params);
}