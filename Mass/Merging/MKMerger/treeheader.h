typedef struct treeparameters{
	int myid,nid;
	float omep,omeplam,omepb;
	float hubble,npow;
	float boxsize,bias;
	float fNL,gNL;
	float amax,anow;
	double pmass;
	float zinit,astep,rsmooth,rth;
	float theta,sphere_radius;
	int nstep,iseed;
	int nx,ny,nz;
} treeparameters;
