#ifdef DEFINE_COMMON_SIMPARAM
#define EXTERN
#else
#define EXTERN extern
#endif
EXTERN int nmesh; 
EXTERN int nxmesh,nymesh,nzmesh; 
EXTERN int nspace0;

EXTERN float boxsize_hinv_mpc;
EXTERN float npow0; 
EXTERN float omega; 
EXTERN float omegabaryon; 
EXTERN float lambda; 
EXTERN float hubble_parameter; 
EXTERN float bias_8mpc; 
EXTERN float rsmooth; 
EXTERN float rtophat_mesh;
EXTERN int powreadflag;
EXTERN char powfilename[190];
EXTERN char inpapkfilename[190];

EXTERN float z_init; 
EXTERN float delta_a;
EXTERN float anow;
EXTERN float amax;
EXTERN int nstep; 
EXTERN int stepnum; 
EXTERN int nskip;

EXTERN float theta;
EXTERN int L, L2; 
EXTERN float SPHERE_RADIUS, particle_radius;

EXTERN float ken;
EXTERN float ktot;
EXTERN float const0;
EXTERN float poten0;
EXTERN float rth;

EXTERN int seed;

EXTERN char rvfilename[80], rvprefix[80];


EXTERN char camerafilename[80];

#undef EXTERN
