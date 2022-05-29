#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "rwalk.h"


int write_data(Sys_Params *params, int *grid, int *visited_list, double weigth, FILE *stream) {
    int N = params->N;
    int L = params->L;

    // Endpoints of walk.
    int end_X, end_Y;
    end_X = visited_list[N-1] % L;
    end_Y = (int)(visited_list[N-1] / L);

    // Walk starts at (N,N). 
    int S2 = (N-end_X)*(N-end_X) + (N-end_Y)*(N-end_Y);

    fprintf(stream, "%d %8e\n", S2, weigth);

    return 1;
}

