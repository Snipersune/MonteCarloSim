#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <math.h>
#include "ising.h"

#define NVAR 5

void sum_scorr(Par *par, FILE *stream){
  
  int if_eof = 0;
  int nbl = 0;

  int L = par->L;
  double *s_corr = malloc(L*sizeof(double));
  double *temp = malloc(L*sizeof(double));

  // Reset s_corr
  for (int i = 0; i < L; i++){
    s_corr[i] = 0.0;
  }
  
  // Run until EOF is reached
  while (1) {
    for (int i = 0; i < L; i++) {
      if (fscanf(stream, "%lf", &temp[i]) != EOF) {
        s_corr[i] += temp[i];
      }
      else {
        if_eof = 1;
        break;
      }
    }

    fscanf(stream, "\n");

    // Check whether EOF has been reached.
    if (if_eof) {
      break;
    }
    nbl++;
  }

  // Output mean of correlation func
  fprintf(stdout, "%d %5e ", par->L, par->t);
  for (int i = 0; i < L; i++) {
    fprintf(stdout, "%8e ", s_corr[i] /= nbl);
  }
  fprintf(stdout, "\n");

  free(temp);
  free(s_corr);
}

int fread_quants(Par *par, double *vsum, double *v2sum, FILE *stream)
{
  int i, nbl = 0;
  double v[NVAR];

  for (i = 0; i < NVAR; i++)
    vsum[i] = v2sum[i] = 0.0;

  // This has to be changed when reading formatted data files.
  while (fscanf(stream, "%lf %lf %lf %lf %lf\n", &v[0], &v[1], &v[2], &v[3], &v[4]) != EOF) {
    for (i = 0; i < NVAR; i++) {
      vsum[i] += v[i];
      v2sum[i] += v[i]*v[i];
    }
    nbl++;
  }
  return nbl;
}

void result(Par *par, double *vsum, double *v2sum, int nblock)
{
  double T2 = par->t*par->t;
  double L2 = par->L*par->L;

  double E, E2, M, M2, M4;

  E = vsum[0]/nblock; E2 = vsum[1]/nblock;
  M = vsum[2]/nblock; M2 = vsum[3]/nblock; M4 = vsum[4]/nblock;

  double M_2 = v2sum[2]/nblock;

  // Heat capacity per spin
  double c = 1.0/(T2)*(E2 - E*E);
  c /= nblock*L2;

  // Standard deviation in estimated mean of m 
  double m_err = sqrt(1.0/(nblock-1)*(M_2 - M*M))/L2;

  // Binder's cumulant
  double Q = M2*M2/M4;

  fprintf(stdout, "%d %8f %8e %8e %8e %8e\n", par->L, par->t, c, M/L2, m_err, Q);
}


int main(int argc, char *argv[])
{
  int i, nblock;
  char fname[FNAMESIZE];

  double vsum[NVAR], v2sum[NVAR];
  double *scorr;

  Par *par = malloc(sizeof(Par));

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
    fscanf(stream, "%d %lf %d %d %d %d\n", &par->L, &par->t, &par->ntherm, &par->nblock, &par->nsamp, &par->seed);

#ifndef SCORR
  nblock = fread_quants(par, vsum, v2sum, stream);
  if (nblock)
    result(par, vsum, v2sum, nblock);
#else
  sum_scorr(par, stream);
#endif

  fclose(stream);
  }

  free(par);
}