
#ifdef SSFAST
#define SS
#endif

#define FNAMESIZE 64


typedef struct Par {
  int L;
  double t;
  int ntherm, nblock, nsamp, seed;
} Par;

typedef struct PhysProps {
  double E;
  double E2;
  double M;
  double M2;
  double M4;
} PhysProps;

// In config.c
extern int write_config(Par *par, int *spin, char *fname);
extern int read_config(Par *par, int *spin, char *fname);

extern int write_data(Par *par, PhysProps *q, char *fname);
extern int read_data(Par *par, PhysProps *q, char *fname);

