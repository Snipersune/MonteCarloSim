#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <time.h>
#include <math.h>
#include "ran.h"
#include <g2.h>
#include <g2_PS.h>
#include <g2_X11.h>
#include <time.h>

#include "ising.h"
#include "my_queue.h"

# define PLUS(x, L) ((x) == (L) - 1 ? 0 : x + 1)
# define MINUS(x, L) ((x) == 0 ? (L) - 1 : x - 1)

#define kb 1.0              // 1.380649e-23
#define J 1.0

// Global variables
char fname[FNAMESIZE];	// Name for config files
double exp_vals[5];     // To store possible exponential values in
double p;              // Cluster acceptance probability

#ifdef PLOT
#define PLOT_VALUE 1
#else
#define PLOT_VALUE 0
#endif
int if_plot = PLOT_VALUE;


void timeDelay(double t_in_sec){
  clock_t t = clock();
  while(clock() < t + (long)(t_in_sec*(double)CLOCKS_PER_SEC)){
    ;
  }
}

void init_tables(Par *par)
{
  extern double exp_vals[5];
  for (int i = 0; i < 5; i++){
    exp_vals[i] = exp(-1.0/(kb*par->t)*J*(4.0*i - 8.0));
  } 
}

void PhysProps_reset(PhysProps * q){
  q->magn = 0;
  q->E = 0;
  q->E2 = 0;
}

double getLocalE(int *spin, int i, int j, int L){

  double E;
  int s_i = spin[i*L + j];

  // Left right up down neigbours
  int nbors_sum = spin[i*L + MINUS(j,L)] + spin[i*L + PLUS(j,L)]
                  + spin[MINUS(i,L)*L + j] + spin[PLUS(i,L)*L + j];
                   
  E = -1.0*J*nbors_sum*s_i;

  return E;
}

double getDeltaE(int *spin, int i, int j, int L){
  double E_orig, E_flipped;

  E_orig = getLocalE(spin, i, j, L);
  E_flipped -= E_orig;

  return E_flipped - E_orig;
}

double getAcceptProb(double t, double dE){

  double alpha = exp(-1.0/(kb*t)*dE);

  if (alpha > 1.0)
    return 1.0;

  return alpha;
}

double getAcceptProbFast(int *spin, int i, int j, int L){

  double alpha;
  extern double exp_vals[5];

  int s_i = spin[i*L + j];

  // Left right up down neigbours
  int nbors_sum = spin[i*L + MINUS(j,L)] + spin[i*L + PLUS(j,L)]
                  + spin[MINUS(i,L)*L + j] + spin[PLUS(i,L)*L + j];

  // Can be -4, -2, 0, 2, 4                 
  nbors_sum = nbors_sum*s_i;

  int table_index = 2 + nbors_sum/2;
  alpha = exp_vals[table_index];

  if (alpha > 1.0)
    return 1.0;

  return alpha;
}

void grow_cluster(Par *par, int *spin, int *nbor_loc, int loc, int state){

  int i = loc/par->L;
  int j = loc%par->L;

  // Check neighbours left right up down
  nbor_loc[0] = i*par->L + MINUS(j,par->L);
  nbor_loc[1] = i*par->L + PLUS(j,par->L);
  nbor_loc[2] = MINUS(i,par->L)*par->L + j;
  nbor_loc[3] = PLUS(i,par->L)*par->L + j;

  int nbor;
  for (int i = 0; i < 4; i++){
    nbor = nbor_loc[i];
    if (spin[nbor] == state && dran() < p){
      insert(nbor);
      spin[nbor] *= -1;
    }
  }
}

double measure(Par *par, PhysProps *q, int *spin){
  double E = 0;
  double magn = 0;

  int L = par->L;
  for (int i = 0; i < L; i++){
    for (int j = 0; j < L; j++){
      
      // Add magnetization of current spin.
      magn += spin[i*L + j];

      // Add local energy contribution. 
      E += getLocalE(spin, i, j, L);
    }
  }

  q->magn += fabs(magn);
  q->E += E/2.0;
  q->E2 += E*E/4.0;

  return E/2.0;
}

void result(Par *par, PhysProps *q, int ntot, int if_final){
  int L2 = par->L*par->L;
  double t = par->t;
  double magn_mean, e_mean, e2_mean, cv;

  magn_mean = q->magn/(double)ntot;
  
  e_mean = q->E/(double)ntot;
  e2_mean = q->E2/(double)ntot;

  cv = 1.0/(kb*t*t)*(e2_mean - e_mean*e_mean);

  magn_mean /= (double)(L2);
  e_mean /= (double)(L2);
  e2_mean /= (double)(L2);
  cv /= (double)(L2);

  if (if_final){
    write_data(par, magn_mean, e_mean, cv, fname);
    printf("  --------  --------  --------\n");
  }
  printf(" %8f  %8f  %8f \n", e_mean, cv, magn_mean);
}

void displaySpin(Par *par, int *spin, int id, int size){
  static int nCalls = 0;

  char plotname[FNAMESIZE] = "plots/spin_";
  char n[10];
  sprintf(n, "%3.3d", nCalls);
  strcat(n, ".ps");
  strcat(plotname, n);

  if (nCalls > 50){
    return;
  }
  int L = par->L;

  int id2;
 
  id2 = g2_open_EPSF_CLIP(plotname, size, size);

  double scale = (double)size/L; 
  
  g2_set_line_width(id, 1);
  g2_set_line_width(id2, 1);
  
  int pen_lg = g2_ink(id, 0.8,0.8,0.8);
  int pen_lg2 = g2_ink(id2, 0.8,0.8,0.8);
  for (int i = 0; i < L; i++){
    for (int j = 0; j < L; j++){

      if (spin[i*L + j] == 1){
        g2_pen(id, pen_lg);
        g2_pen(id2, pen_lg2);
      }
      else{
        g2_pen(id, 1);
        g2_pen(id2, 1);
      }
      g2_filled_rectangle(id, scale*i, scale*j, scale*(i+1), scale*(j+1));
      g2_filled_rectangle(id2, scale*i, scale*j, scale*(i+1), scale*(j+1));    
    }
  }
  g2_close(id2);

  timeDelay(0.02);

  nCalls++;
}

void writeTimeCorrelationFunc(Par *par, double *E_vals, double E_mean){
  
  char filename[32];
  sprintf(filename, "data/corr_%d_%.3f", par->L, par->t);
  FILE *fp = fopen(filename, "w");
  if (fp == NULL){
    printf("Could not open correlation time file!\n");
    return;
  }
  fprintf(fp, "t  C_e\n");
  
  double L2 = (double)par->L*par->L;
  //double e_mean = E_mean/L2;
  //double e2 = e_mean*e_mean;

  double E2 = E_mean*E_mean;

  double E0Et_mean, E0Et;
  int n_vals = par->nblock*par->nsamp;

  int t0;
  for (int t = 0; t < 201; t++){
    t0 = 0;
    E0Et = 0;
    while(t0 + t < n_vals){
      E0Et += E_vals[t0]*E_vals[t0 + t];
      t0++;
    }
    E0Et_mean = E0Et/((double)t0);
    fprintf(fp, "%d %f\n", t, (E0Et_mean - E2)/((double)L2));
  }
  fclose(fp);
}

#ifndef CLU
int update(Par *par, int *spin) {

  static int if_called = 0;
  int accept = 0;
  
  double alpha, t = par->t;
  int L = par->L, L2 = L*L;

  int x11_id, size;
  size = 600;
  if (!if_called && if_plot){
    x11_id = g2_open_X11(size, size);
  }

  double r;
  for (int i = 0; i < L; i++){
    for (int j = 0; j < L; j++){

      if(!if_called && if_plot){
        displaySpin(par, spin, x11_id, size);
      }

      // Get probability of accepting state change.
      alpha = getAcceptProbFast(spin, i, j, L);

      // If change state.
      r = dran();
      if (r <= alpha){
        spin[i*L + j] *= -1.0;
        accept++;
      }
    }
  }

  if (!if_called && if_plot){
    g2_close(x11_id);
  }
  if_called++;

  return accept;
}
#endif

#ifdef CLU
#ifdef TRI
#define NNN 6
#else
#define NNN 4
#endif

int update(Par *par, int *spin){

  static int if_called = 0;

  int x11_id, size;
  size = 600;
  if (!if_called && if_plot){
    x11_id = g2_open_X11(size, size);
  }

  int accept = 0;
  // Get random start location
  int start_loc = (int)(dran()*par->L*par->L);
  int start_state = spin[start_loc];

  // Insert in queue to initialize loop
  insert(start_loc);
  spin[start_loc] *= -1;

  int nbor_loc[4];
  int new_loc;
  while(!isEmpty()){
    if(!if_called && if_plot) {
        displaySpin(par, spin, x11_id, size);
      }

    new_loc = removeData();
    accept++;
    grow_cluster(par, spin, nbor_loc, new_loc, start_state);
  }
  if (!if_called && if_plot){
    g2_close(x11_id);
  }
  if_called++;
  
  return accept;
}
#endif


int mc(Par *par, int *spin)
{
  int i, iblock, isamp, istep, ntherm = par->ntherm;
  double t = par->t, acc, accept = 0.0, L2 = par->L * par->L;

  double *E_vals;

#ifdef CORR
  E_vals = malloc((par->nblock*par->nsamp)*sizeof(double)); // To store E(t) for correlation time
#endif

  PhysProps q_tot = {0, 0, 0}, q = {0, 0, 0};

  // *** Read in the configuration for the present parameters if already present.
  if (read_config(par, spin, fname))
    ntherm = 0;

  // *** Write out information about the run: size, temperature,...
#ifdef CLU
  printf("2D Ising model with the Wolff cluster algorithm.\n");
#else
  printf("2D Ising model with the Metropolis algorithm.\n");
#endif

  printf("\n====    %d x %d     T = %g    ====\n", par->L, par->L, par->t);
  printf("\nntherm  nblock   nsamp   seed\n");
  printf(" %5d   %5d   %5d   %d\n", ntherm, par->nblock, par->nsamp, par->seed);

  printf("\n energy      cv        magn     \n");

#ifdef CLU
  SetQueueMAXSize(L2);
  if (!init_queue()) {
    return 0;
  }
  p = 1.0 - exp(-2.0/(kb*par->t));
#endif

  // Thermalize the system 
  for (i = 0; i < ntherm; i++)
    update(par, spin);

  double E;
  for (iblock = 0; iblock < par->nblock; iblock++) {
    for (isamp = 0; isamp < par->nsamp; isamp++) {

      accept += update(par, spin);
      E = measure(par, &q, spin);

#ifdef CORR
  E_vals[iblock*par->nsamp + isamp] = E;
#endif
    }
    
    // Fix this (5): write_config(par, spin, fname);
    // write_config(par, spin, fname);

    result(par, &q, par->nsamp, 0);

    q_tot.magn += q.magn;
    q_tot.E += q.E;
    q_tot.E2 += q.E2;
    PhysProps_reset(&q);
  }

  result(par, &q_tot, par->nblock*par->nsamp, 1);
  write_config(par, spin, fname);

  acc = accept * 100.0 / (L2 * par->nblock * par->nsamp);
  printf("\nAcceptance: %5.2f\n", acc);

#ifdef CORR
  double E_mean = q_tot.E/((double)par->nblock*par->nsamp);
  writeTimeCorrelationFunc(par, E_vals, E_mean);
  free(E_vals);
#endif

#ifdef CLU
  free_queue();
#endif

  return 1;
}


int initialize_mc(Par *par, int *spin) {
  int i, L2 = par->L * par->L;

  if (!par->L) {
    printf("Give system size N!\n");
    return 0;
  }


  init_ran(par->seed);

  double exp_vals[5];
  init_tables(par);

  sprintf(fname, "%3.3d_%5.3f", par->L, par->t);

  return mc(par, spin);
}


int read_args(Par *par, char *arg)
{
  int st;
  static int *spin = NULL;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_mc(par, spin);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "L")) {
    int i, L2;
    par->L = strtod(s, NULL);



    L2 = par->L * par->L;
    spin = realloc(spin, L2 * sizeof(int));
    for (i = 0; i < L2; i++)
      spin[i] = 1;
#ifdef CLU
    par->ntherm = 1000;
#else
    par->ntherm = L2;
#endif

    return 1;
  }

  if (!strcmp(arg, "T")) {
    par->t = strtod(s, NULL);
    return 1;
  }


  if (!strcmp(arg, "ntherm")) {
    par->ntherm = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "nblock")) {
    par->nblock = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "nsamp")) {
    par->nsamp = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "seed")) {
    par->seed = strtol(s, NULL, 0);
    return 1;
  }

  if (!strcmp(arg, "read")) {
    return read_config(par, spin, s);
  }


  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}


int main(int argc, char *argv[])
{
  int i, iarg;
  Par *par = malloc(sizeof(Par));

  par->L = 0;
  par->t = 2.26;
  par->nblock = 1;
#ifdef CLU
  par->nsamp = 1000;
#else
  par->nsamp = 10000;
#endif

  if (argc == 1) {
    printf("Usage: %s L=16 T=2.26\n", argv[0]);
    printf("Optional arguments (with defaults) nblock=%d nsamp=%d ntherm=%d seed=%d\n",
	   par->nblock, par->nsamp, par->ntherm, par->seed);
    exit(EXIT_SUCCESS);
  }

  // Interpret the commands given in the argument list.
  for (iarg = 1; iarg < argc; iarg++)
    if (!read_args(par, argv[iarg]))
      exit(EXIT_FAILURE);

  free(par);
  return 0;
}
