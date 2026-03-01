void totalken_(float *,float *,float *,int *,float *,float *,float *,int *);
void readpar_(int *, float *,float *,float *);
void pmseed_(int *,pmparticletype *, int *, int *,
		float *, float *, float *,float *,float *,float *,float *,float *,float *);
void subpmseed_(long *, pmparticletype *,int *,int *,int *,
		float *,int *,float *,float *,float *,float *, float *,float *,float *,
		int *,float *,float *);
void setinnerparticles_(long *, pmparticletype *,int *,int *,int *,
		float *,int *,float *,float *,float *,float *, float *,float *,float *,
		int *,float *,float *);
void setinnerparticlesforce_(long *, pmparticletype *,int *,int *,int *,
		float *,int *,float *,float *,float *, float *,float *,float *,
		int *,float *,float *);
void getinnerparticlesforce_(long *, pmvrparticletype *,int *,int *,int *,
		float *,int *,float *,float *,float *, float *,float *,float *,
		int *,float *,float *,float *,float *,float *);
void setboundaryparticles_(long *, pmparticletype *,int *,int *,int *,
		float *,int *,float *,float *,float *,float *, float *,float *,float *,
		int *,float *,float *);
void setboundaryparticlesforce_(long *, pmparticletype *,int *,int *,int *,
		float *,int *,float *,float *,float *,float *, float *,float *,
		int *,float *,float *,float *,float *,float *);
#ifdef PMSEEDFORCE
void pmseedforce_(long *, pmparticletype *, int *, int *,
		float *, float *, float *,float *,float *,int *,float *,
		float *, float *,float *,float *, float *, float *,float *,
		float * );
#else
void pmseedforce_(long *, pmparticletype *, int *, int *,
		float *, float *, float *,float *,float *,float *,float *,float *, float *);
#endif
void fda4sup_(int *,int *,int *,int *);
void npdata_(int *,int *,int *,int *);
void fftwinit_(int *,int*,int*,int*,int*);

void rdatafromdisk_(int *,int *,int *,int *,float *,float *,float *,
		float *,float *,float *,float *,float *,float *,float *,float *,
		float *,float *,float *,float *,float *,int *,int *,int *,int *,
		int *,int *);
void fftforward_(float *);
void fftbackward_(float *);
void tsc_(float *,int *,int *,pmparticletype *, int *,
		int *,int *,float *,int *, int *);
void psolver_(float *,float *,float *);
/*
void fda4gf_(pmparticletype *, int *,
		int *,float *,float *,float *,float *,int *,int *,int *,int *,
		float *,float *,float *,float *,float *,
		float *,float *,int *,int *,int *,int *,int *,int *,int *,int *,int *,
		float *,float *);
		*/
void fda4inner_(pmparticletype *,int *, int *,
		int *,float *,int *,float *,float *,float *,
		float *,float *,float *,float *,float *,
		int *,float*,float*,int *,int *);
void fda4boundary_(pmparticletype *,int *, int *,
		int *,float *,int *,float *,float *,float *,
		float *,float *,float *,float *,float *,
		int *,float*,float*,int *,int *, int *);
void fda4_(pmparticletype *,int *, int *,
		int *,float *,int *,float *,float *,float *,
		float *,float *,float *,float *,float *,
		int *,float*,float*,int *,int *);
void wdata2disk_(int *,int *,int *, int *,float *,float *,float *,
		float *,float *,float *,float *,float *,float *,float *,float *,float *,
		float *,float *,float *,int *,int *,int *,int *,int *);
