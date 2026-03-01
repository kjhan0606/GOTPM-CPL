void write_head(FILE *, SimParameters);
SimParameters read_head(FILE*);
#define MAX_LINE_LENGTH 200
#define S_FLOAT   "g"
#define S_DOUBLE   "lf"
#define S_INT     "d"
#define S_LONG     "ld"
#define S_STRING  "s"
#define SET     "define "

#define P_Starting  "#Start of the Ascii Header of the GalMTree\n"
#define P_IC        SET"INITIAL CONDITION   = %"S_INT"\n"
#define P_Nid       SET"Nid                 = %"S_INT"\n"
#define P_Myid      SET"Myid                = %"S_INT"\n"
#define P_Omep      SET"OmegaMatter0        = %"S_FLOAT"\n"
#define P_Omepb     SET"OmegaBaryon0        = %"S_FLOAT"\n"
#define P_Omeplamb  SET"OmegaLambda0        = %"S_FLOAT"\n"
#define P_Omei      SET"OmegaMatterI        = %"S_FLOAT"\n"
#define P_fNL       SET"fNL                 = %"S_FLOAT"\n"
#define P_gNL       SET"gNL                 = %"S_FLOAT"\n"
#define P_Hubble    SET"Hubble              = %"S_FLOAT"\n"
#define P_nPS       SET"nPS                 = %"S_FLOAT"\n"
#define P_Boxsize   SET"Boxsize(Mpc/h)      = %"S_FLOAT"\n"
#define P_Bias      SET"Bias                = %"S_FLOAT"\n"
#define P_Amax      SET"Amax                = %"S_FLOAT"\n"
#define P_Anow      SET"Anow                = %"S_FLOAT"\n"
#define P_Nx        SET"Nx                  = %"S_INT"\n"
#define P_Ny        SET"Ny                  = %"S_INT"\n"
#define P_Nz        SET"Nz                  = %"S_INT"\n"
#define P_Zinit     SET"Z_init              = %"S_FLOAT"\n"
#define P_Astep     SET"Astep               = %"S_FLOAT"\n"
#define P_Nstep     SET"Nstep               = %"S_INT"\n"
#define P_Iseed     SET"Iseed               = %"S_INT"\n"
#define P_XYZSHIFT  SET"XYZ_Shift_Flag      = %"S_INT"\n"
#define P_PTSIZE    SET"Particle_TYPE_Size  = %"S_INT"\n"
#define P_Closing   "#End of Ascii Header\n"
#define P_NULL      "##################################################\n"

#define MPI_DATA(frw,wp,sp,simpar){\
 	ncnt += frw(wp,P_Myid,sp simpar.myid);\
	ncnt += frw(wp,P_Nid,sp simpar.nid);\
}
#define CORE(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### INITIAL CONDITION: 1 (ZA)/ 2 (2'nd order LPT)\n");\
 	ncnt += frw(wp,P_IC,sp simpar.IC);\
	ncnt += frw(wp,P_Omep,sp simpar.omep);\
	ncnt += frw(wp,P_Omepb,sp simpar.omepb);\
	ncnt += frw(wp,P_Omeplamb,sp simpar.omeplam);\
	ncnt += frw(wp,P_fNL,sp simpar.fNL);\
	ncnt += frw(wp,P_gNL,sp simpar.gNL);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### Hubble parameter is in 100km/sec/Mpc.\n");\
	ncnt += frw(wp,"### nPS is the power spectral index.\n");\
	ncnt += frw(wp,"### Boxsize should be in unit of Mpc/h.\n");\
	ncnt += frw(wp,"### Bias factor is inverse of sigma_8.\n");\
	ncnt += frw(wp,P_Hubble,sp simpar.hubble);\
	ncnt += frw(wp,P_nPS,sp simpar.npow);\
	ncnt += frw(wp,P_Boxsize,sp simpar.boxsize);\
	ncnt += frw(wp,P_Bias,sp simpar.bias);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Amax,sp simpar.amax);\
	ncnt += frw(wp,P_Anow,sp simpar.anow);\
	ncnt += frw(wp,P_Astep,sp simpar.astep);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### Nx, Ny, and Nz should be same & nspace should be int.\n");\
	ncnt += frw(wp,P_Nx,sp simpar.nx);\
	ncnt += frw(wp,P_Ny,sp simpar.ny);\
	ncnt += frw(wp,P_Nz,sp simpar.nz);\
	ncnt += frw(wp,P_Nstep,sp simpar.nstep);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Iseed,sp simpar.iseed);\
	ncnt += frw(wp,P_NULL);\
}
#define DEFAULT_PARAMS(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_Starting);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_Closing);\
}

#define PARAMS(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_Starting);\
	CORE(frw,wp,sp,simpar);\
	simpar.ptypesize = sizeof(pmparticletype);\
	ncnt += frw(wp,P_PTSIZE,sp simpar.ptypesize);\
	ncnt += frw(wp,P_Closing);\
}

#define FILE_HEADER(frw,wp,sp,simpar) {\
	ncnt += frw(wp,P_Starting);\
	ncnt += frw(wp,P_NULL);\
	MPI_DATA(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_NULL);\
	CORE(frw,wp,sp,simpar);\
	simpar.ptypesize = sizeof(pmparticletype);\
	ncnt += frw(wp,P_Np,sp simpar.np);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Closing);\
}
