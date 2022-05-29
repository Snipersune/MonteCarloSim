#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ran.h"
#include "perc.h"
#include "my_queue.h"


#ifndef THREE_D
#define DIM 2
#else
#define DIM 3
#endif

#define PLUS(x, L) ((x) == (L) - 1 ? 0 : x + 1)
#define MINUS(x, L) ((x) == 0 ? (L) - 1 : x - 1)

#define VISITED 10 
#define XMARK 100
#define YMARK 1000
#define ZMARK 10000

// Global variables
char fname[FNAMESIZE];	// Name for config files
int D = DIM;


// Relative position of nearest neighbours of a site. 
Point step[6] = {
    {0, 1, 0},         // Up
    {0, -1, 0},        // Down
    {-1, 0, 0},        // Left
    {1, 0, 0},         // Right
    {0, 0, 1},         // In
    {0, 0, -1},        // Out
};


void initialize_grid(Sys_Params *params, int *grid) {
    int Lx = params->Lx;
    int Ly = params->Ly;
    int Lz = params->Lz;

    int L2 = Lx*Ly;
    
    double p = params->p;

    int pos_index;
    double r;
    for (int z = 0; z < Lz; z++) {
        for (int y = 0; y < Ly; y++) {
            for (int x = 0; x < Lx; x++) {
                pos_index = z*L2 + y*Ly + x;
                r = dran();
                if (r < p) {
                    grid[pos_index] = 1;
                } 
                else {
                    grid[pos_index] = 0;
                }
            }
        }
    }
}

void getNNbors(int Lx, int Ly, int Lz, int x, int y, int z, int *nbor) {
    int Lxy = Lx*Ly;

    // Get neighbours up down left right
    nbor[0] = Lxy*z + Ly*PLUS(y,Ly) + x;
    nbor[1] = Lxy*z + Ly*MINUS(y,Ly) + x;
    nbor[2] = Lxy*z + Ly*y + MINUS(x,Lx);
    nbor[3] = Lxy*z + Ly*y + PLUS(x,Lx);

    // If 3D get behind infront also
#ifdef THREE_D
    nbor[4] = Lxy*PLUS(z,Lz) + Ly*y + x;
    nbor[5] = Lxy*MINUS(z,Lz) + Ly*y + x;
#endif

}

// Returns 1 if coord > L, -1 if coord < 0, else 0;
int region_change(int coord, int L) {
    if (coord < 0)
        return -1;
    else if (coord > L-1)
        return 1;
    
    return 0;
}

// Returns marker change depending on region change
int get_marker_change(int Lx, int Ly, int Lz, int nx, int ny, int nz) {

    // To store region change
    int rc;

    // Get marker change
    if ((rc = region_change(nx, Lx)) != 0)        // Check x-axis
        return rc*XMARK;

    else if ((rc = region_change(ny, Ly)) != 0)   // y-axis   
        return rc*YMARK;

    else if ((rc = region_change(nz, Lz)) != 0)   // z-axis
        return rc*ZMARK;
    
    return 0;
}

int check_percolation(Sys_Params *params, int *grid, int start_x, int start_y, int start_z) {
    int Lx, Ly, Lz, Lxy; 
    
    Lx = params->Lx;
    Ly = params->Ly;
    Lz = params->Lz;
    Lxy = Lx*Ly;

    int start_index;
    start_index = Lxy*start_z + Ly*start_y + start_x;

    // Set visited and instert into queue
    grid[start_index] = VISITED;
    insert(start_index);

    int nbors[2*DIM];

    int curr_index;
    int x, y, z;
    while (!isEmpty()) {

        // Get site from queue
        curr_index = removeData();

        // Get corresponding x, y, z coordinate
        z = curr_index / Lxy;
        y = (curr_index-Lxy*z) / Ly;
        x = (curr_index-Lxy*z) % Ly;

        // Get neighbours
        getNNbors(Lx, Ly, Lz, x, y, z, nbors);

        // Go through neighbours
        for (int i = 0; i < 2*DIM; i++) {

            // If occupied
            if (grid[nbors[i]] != 0) {
                int nx, ny, nz;

                // Get neighbour position
                nx = x + step[i].X;
                ny = y + step[i].Y;
                nz = z + step[i].Z;

                int marker_change;
                
                // Store marker change for neighbour if region change
                marker_change = get_marker_change(Lx, Ly, Lz, nx, ny, nz);
                
                // If unvisited
                if (grid[nbors[i]] == 1) {

                    grid[nbors[i]] = grid[curr_index];     // Inherit marker
                    insert(nbors[i]);                       // Add to queue

                    // Apply marker change
                    grid[nbors[i]] += marker_change;
                }

                // If previously visited
                else {
                    // Store marker assigned when visited
                    int test_marker;
                    int assigned_marker = grid[nbors[i]];

                    
                    // Apply marker change
                    test_marker = grid[curr_index] + marker_change;

                    // If not same marker, site is reached in different region and perculation is found
                    if (test_marker != assigned_marker) {
                        return 1;
                    }
                }
            }                     
        }
    }

    // No percolation found
    return 0;
}

int if_percolation(Sys_Params *params, int *grid) {

    int Lx = params->Lx;
    int Ly = params->Ly;
    int Lz = params->Lz;
    int Lxy = Lx*Ly;

    // Check xz row or plane (3D) for y=0.
    for (int z = 0; z < Lz; z++) {
        for (int x = 0; x < Lx; x++) {

            // If site is occupied but not visited, start percolation check.
            if (grid[Lxy*z + x] == 1) {
                if (check_percolation(params, grid, x, 0, z)) {
                    clearQueue();
                    return 1;
                }
            }   
        }
    }

    // Check yz row or plane (3D) for x=0.
    for (int z = 0; z < Lz; z++) {
        for (int y = 0; y < Ly; y++) {

            // If site is occupied but not visited, start percolation check.
            if (grid[Lxy*z + Ly*y] == 1) {
                if (check_percolation(params, grid, 0, y, z)) {
                    clearQueue();
                    return 1;
                }
            }   
        }
    }

#ifdef THREE_D
    // Check xy plane for z=0.
    for (int y = 0; y < Ly; y++) {
        for (int x = 0; x < Lx; x++) {

            // If site is occupied but not visited, start percolation check.
            if (grid[Ly*y + x] == 1) {
                if (check_percolation(params, grid, x, y, 0)) {
                    clearQueue();
                    return 1;
                }
            }   
        }
    }
#endif

    return 0;
}

int perc(Sys_Params *params, char* fp) {
    
    int N = params->N;

    int Lx, Ly, Lz, Lxy;
    Lx = params->Lx;
    Ly = params->Ly;
    Lz = params->Lz;
    Lxy = Lx*Ly;

    double p;
    p = params->p;

    fprintf(stdout, "Check percolation for system of size (Lx,Ly,Lz)=(%d,%d,%d) and p=%4lf.\n", Lx, Ly, Lz, p);
    fprintf(stdout, "Simulating for %d systems.\n", N);

    int *grid = malloc(Lz*Lxy*sizeof(int));

    if (grid == NULL)
        return 0;

    // Initialize queue used when checking percolation
    setQueueMAXSize(Lx*Ly*Lz);
    init_queue();

    // Check perculation for N systems.
    int n_percs = 0;
    for (int i = 0; i < N; i++) {

        initialize_grid(params, grid);
        n_percs += if_percolation(params, grid);
    }

    write_data(params, n_percs, fname);

    free_queue();
    free(grid);

    fprintf(stdout, "Simulation complete!\n\n");
    return 1;
}

int initialize_perc(Sys_Params *params) {

    if (!params->Lx) {
        printf("Give system size L!\n");
        return 0;
    }
    
    init_ran(params->seed);

    sprintf(fname, "%3.3d_%.3f", params->Lx, params->p);

    char path[FNAMESIZE + 5] = "data/";

    strcat(path, fname);
    FILE *fp;
    // If file already exists, append to existing data. Otherwise write parameter data.
    if (fp = fopen(path, "r")) {
        fclose(fp);
    }
    else if (fp = fopen(path, "w")) {
        fprintf(fp, "Lx=%d Ly=%d Lz=%d p=%.8e\n", params->Lx, params->Ly, params->Lz, params->p);
        fclose(fp);
    }
    else {
        fprintf(stdout, "Could not open file to write data. Program terminated!\n");
        return 0;
    }

    return perc(params, fname);
}


int read_args(Sys_Params *params, char *arg) {
  char *ptr;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_perc(params);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "N")) {
    params->N = strtol(s, &ptr, 10);
    if (ptr == arg){
            fprintf(stderr, "Non-valid value of N entered, must be integer valued. Program exit!\n");
            return 0;
        }
    return 1;
  }

  if (!strcmp(arg, "L")) {
    params->Lx = strtol(s, &ptr, 10);
    params->Ly = params->Lx;
    params->Lz = 1;
#ifdef THREE_D
    params->Lz = params->Lx;
#endif
    if (ptr == arg){
            fprintf(stderr, "Non-valid value of L entered, must be integer valued. Program exit!\n");
            return 0;
        }
    return 1;
  }

  if (!strcmp(arg, "p")) {
    params->p = strtod(s, &ptr);
    if (ptr == arg || params->p < 0){
            fprintf(stderr, "Non-valid value of p entered, must be a positive valued real number. Program exit!\n");
            return 0;
        }
    return 1;
  }

  if (!strcmp(arg, "seed")) {
    params->seed = strtol(s, &ptr, 10);
    if (ptr == arg){
            fprintf(stderr, "Non-valid value of W entered, must be integer valued. Program exit!\n");
            return 0;
        }
    return 1;
  }

  fprintf(stderr, "No such variable name: '%s'.\n", arg);
  return 0;
}


int main(int argc, char* argv[]){

    Sys_Params *params = malloc(sizeof(Sys_Params));

    params->N = 0;
    params->Lx = 0;
    params->Ly = 0;
    params->Lz = 0;
    params->p = 0.0;

    // If no command line arguments
    if (argc == 1) {
        printf("Must enter system size 'L=<VAL>' and number of configurations 'N=<VAL>' for program to run.\n");
        exit(EXIT_SUCCESS);
    }

    // Interpret the commands given in the argument list.
    for (int arg = 1; arg < argc; arg++)
        if (!read_args(params, argv[arg]))
            exit(EXIT_FAILURE);

  free(params);
  return 0;
}