#define FNAMESIZE 64

typedef struct Sys_Params {
    int Lx;         // Size in X direction
    int Ly;         // Size in Y direction   
    int Lz;         // Size in Z direction
    double p;          // Probability of occupancy
    int N;          // Number of systems to be generated.
    int seed;
} Sys_Params;


typedef struct Point {
    int X;
    int Y;
    int Z;
} Point;

int write_data(Sys_Params *params, int n_percs, char *fp);
