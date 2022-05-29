
#define FNAMESIZE 64

typedef struct Sys_Params {
    int L;      // Grid size
    int N;      // Walk length
    int W;      // Number of walks
    int seed;
} Sys_Params;

typedef struct Point {
    int X;
    int Y;
} Point;


int write_data(Sys_Params *params, int *grid, int *visited_list, double weigth, FILE *stream);
