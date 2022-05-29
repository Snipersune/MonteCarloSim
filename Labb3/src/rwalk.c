#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "ran.h"
#include "rwalk.h"

// Global variables
char fname[FNAMESIZE];	// Name for config files

// Contains valid directions depending on last move. 0=UP, 1=DOWN, 2=LEFT, 3=RIGHT
int valid_dirs[4][3] = {
    {0, 2, 3},
    {1, 2, 3},
    {0, 1, 2},
    {0, 1, 3},
};

Point step[4] = {
    {0, 1},         // Up
    {0, -1},        // Down
    {-1, 0},        // Left
    {1, 0}          // Right
};

int get_point_index(Sys_Params *params, Point point) {
    int L = params->L;
    return point.X + point.Y*L;
}

void grid_reset(Sys_Params *params, int *grid, int *visited_list) {
    int N = params->N;

    for (int i = 0; i < N; i++)
        grid[visited_list[i]] = 0;
}

int set_visited(Sys_Params *params, int *grid, Point point) {

    int point_i = get_point_index(params, point);
    grid[point_i] = 1;
    
    return point_i;
}

int if_visited(Sys_Params *params, int *grid, Point point) {

    return grid[get_point_index(params, point)];
}

#ifndef SURR
int perform_rwalk(Sys_Params *params, int *grid, int *visited_list, double *weight) {

    int N = params->N;

    grid_reset(params, grid, visited_list);

    int curr_i;
    Point curr_pos;

    // Start at center of grid. 
    curr_pos.X = N;
    curr_pos.Y = N;
    curr_i = set_visited(params, grid, curr_pos);

    int next_step, prev_step, rand_i;
    
    // Perform first move
    next_step = (int)(4*dran());
    curr_pos.X += step[next_step].X;
    curr_pos.Y += step[next_step].Y;

    // Mark as visited
    curr_i = set_visited(params, grid, curr_pos);
    visited_list[0] = curr_i;

    prev_step = next_step;

    int n_steps = 1;
    for (int i = 1; i < N; i++) {

        // Random next step direction, from current valid options.
        rand_i = (int)(3*dran());

        next_step = valid_dirs[prev_step][rand_i];

        // Perform step.
        curr_pos.X += step[next_step].X;
        curr_pos.Y += step[next_step].Y;
        
        // Check if walk is intersecting.
        if (if_visited(params, grid, curr_pos)) {
            break;
        }

        // Otherwise set visited.
        curr_i = set_visited(params, grid, curr_pos);
        visited_list[i] = curr_i;

        // Update last step.
        prev_step = next_step;
        
        n_steps++;
    }
    *weight = 1.0;

    return n_steps;
}
#else
int perform_rwalk(Sys_Params *params, int *grid, int *visited_list, double *weight) {

    int N = params->N;

    grid_reset(params, grid, visited_list);

    int L = params->L;

    // To store valid directions. Will at max be three different directions.
    int possible_dirs[3];
    int n_dirs;

    int curr_i;
    Point curr_pos;

    // Start at center of grid. 
    curr_pos.X = N;
    curr_pos.Y = N;
    curr_i = set_visited(params, grid, curr_pos);

    int next_step, prev_step, rand_i;    
    
    // Perform first move
    next_step = (int)(4*dran());
    curr_pos.X += step[next_step].X;
    curr_pos.Y += step[next_step].Y;

    // Mark as visited
    curr_i = set_visited(params, grid, curr_pos);
    visited_list[0] = curr_i;

    prev_step = next_step;

    *weight = 1.0;

    Point temp;
    int n_steps = 1;
    for (int i = 1; i < N; i++) {

        // Check valid direction.
        n_dirs = 0;
        int test_dir, test_X, test_Y;
        for (int d = 0; d < 3; d++) {
            
            test_dir = valid_dirs[prev_step][d];

            test_X = step[test_dir].X;
            test_Y = step[test_dir].Y;
            
            temp.X = curr_pos.X + test_X;
            temp.Y = curr_pos.Y + test_Y;

            if (!if_visited(params, grid, temp)) {
                possible_dirs[n_dirs] = test_dir;  
                n_dirs++;
            }
        }

        // Stop if no valid direction is available.
        if (n_dirs == 0) {
            break;
        }

        // Random next step direction, from current valid options.
        rand_i = (int)(n_dirs*dran());
        next_step = possible_dirs[rand_i];

        // Perform step.
        curr_pos.X += step[next_step].X;
        curr_pos.Y += step[next_step].Y;
        curr_i = set_visited(params, grid, curr_pos);
        
        visited_list[i] = curr_i;
        
        // Update last step.
        prev_step = next_step;

        // Update weigth.
        *weight *= (double)n_dirs/3.0;
        
        n_steps++;
    }

    return n_steps;
}

#endif


int rwalk(Sys_Params *params, FILE *stream) {
    int L, N, W;

    L = params->L;
    N = params->N;
    W = params->W;

    fprintf(stdout, "Random walk with N=%d and W=%d.\n", N, W);

    int *grid = calloc(L*L, sizeof(int));
    int *visited_list = calloc(N, sizeof(int));

    int walk_length;
    
    // Make random walk until 'W' walks of length 'N' has been achieved.
    double weigth;
    int num_Nwalks = 0;
    while (num_Nwalks < W){
        walk_length = perform_rwalk(params, grid, visited_list, &weigth);

        // Write data to file and increment number of full walks made.
        if(walk_length == N){
            write_data(params, grid, visited_list, weigth, stream);
            num_Nwalks++;
        }
    }

    free(grid);
    free(visited_list);

    fprintf(stdout, "Simulation complete!\n\n");
    return 1;
}

int initialize_rwalk(Sys_Params *params) {

    if (!params->N) {
        printf("Give walk size N!\n");
        return 0;
    }
    
    init_ran(params->seed);

    sprintf(fname, "%3.3d", params->N);

    char path[FNAMESIZE + 5] = "data/";

    strcat(path, fname);
    FILE *fp;
    // If file already exists, append to existing data. Otherwise write parameter data.
    if(fp = fopen(path, "r")){
        fclose(fp);
        fp = fopen(path, "a");
    }
    else if (fp = fopen(path, "w")) {
        fprintf(fp, "%d %d %d\n", params->N, params->W, params->seed);
    }
    else {
        fprintf(stdout, "Could not open file to write data. Program terminated!\n");
        return 0;
    }

    rwalk(params, fp);
    fclose(fp);

    return 1;
}


int read_args(Sys_Params *params, char *arg) {
  char *ptr;
  char *s;

  if (!strcmp(arg, "run"))
    return initialize_rwalk(params);

  s = strchr(arg, '=');

  if (!s) {
    fprintf(stderr, "Command '%s' not recognized, expected format: '<name>=<value>'.", arg);
    return 0;
  }

  *s++ = '\0';

  if (!strcmp(arg, "N")) {
    params->N = strtol(s, &ptr, 10);
    params->L = 2*params->N + 1;
    if (ptr == arg){
            fprintf(stderr, "Non-valid value of N entered, must be integer valued. Program exit!\n");
            return 0;
        }
    return 1;
  }

  if (!strcmp(arg, "W")) {
    params->W = strtol(s, &ptr, 10);
    if (ptr == arg){
            fprintf(stderr, "Non-valid value of W entered, must be integer valued. Program exit!\n");
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
    params->L = 0;
    params->W = 0;

    // If no command line arguments
    if (argc == 1) {
        printf("Must enter walk length 'N=<VAL>' and number of walks 'W=<VAL>' for program to run.\n");
        exit(EXIT_SUCCESS);
    }

    // Interpret the commands given in the argument list.
    for (int arg = 1; arg < argc; arg++)
        if (!read_args(params, argv[arg]))
            exit(EXIT_FAILURE);

  free(params);
  return 0;
}