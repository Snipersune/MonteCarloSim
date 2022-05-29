#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "rwalk.h"

#define NVAR 2

int read_S2(Sys_Params *params, double *vsum, double *v2sum, FILE *stream) {
    int nbl = 0;
    double v[NVAR];

    for (int i = 0; i < NVAR+1; i++){
        vsum[i] = 0;
        v2sum[i] = 0;
    }

    while (fscanf(stream, "%lf %lf\n", &v[0], &v[1]) != EOF) {
        for (int i = 0; i < NVAR; i++) {
            vsum[i] += v[i];
            v2sum[i] += v[i]*v[i];
        }

        vsum[NVAR] += v[0]*v[1]; 
        v2sum[NVAR] += v[0]*v[0]*v[1]; 
        nbl++;
    }

    return nbl;
}

void result(Sys_Params *params, double *vsum, double *v2sum, int nblock) {
  
    int N = params->N;
    double S2w_mean, S2_mean, S2_squared_mean, S2_var, S2w_err;

    double wsum, w2sum;
    wsum = vsum[1];
    w2sum = v2sum[1];

    // Weighted mean
    S2w_mean = vsum[NVAR]/wsum;

    // Unweighted mean and square mean.
    S2_mean = vsum[0]/nblock;
    S2_squared_mean = v2sum[0]/nblock;


    // Standard deviation in estimated mean of S2
    S2_var = nblock/((double)(nblock-1))*(S2_squared_mean - S2_mean*S2_mean);
    S2w_err = sqrt(S2_var*w2sum/(wsum*wsum));

    fprintf(stdout, "%d %6e %6e\n", N, S2w_mean, S2w_err);
}

int main(int argc, char *argv[]) {

    char fname[FNAMESIZE];
    int nblock;

    double vsum[NVAR+1], v2sum[NVAR+1];
    
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
        fscanf(stream, "%d %d %d\n", &params->N, &params->W, &params->seed);

        nblock = read_S2(params, vsum, v2sum, stream);
        if (nblock)
            result(params, vsum, v2sum, nblock);

        fclose(stream);
    }

  free(params);
}