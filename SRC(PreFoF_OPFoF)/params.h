void write_head(FILE *, SimParameters);
SimParameters read_head(FILE*);
#define MAX_LINE_LENGTH 200
#define S_FLOAT   "g"
#define S_DOUBLE   "lf"
#define S_INT     "d"
#define S_LONG     "ld"
#define S_STRING  "s"
#define SET     "define "

#define P_Starting  "#Start of the Ascii Header of the GOTPM Simulation\n"
#define P_IC        SET"INITIAL CONDITION   = %"S_INT"\n"
#define P_Nid       SET"Nid                 = %"S_INT"\n"
#define P_Myid      SET"Myid                = %"S_INT"\n"
#define P_Omep      SET"OmegaMatter0        = %"S_FLOAT"\n"
#define P_Omepb     SET"OmegaBaryon0        = %"S_FLOAT"\n"
#define P_Omeplamb  SET"OmegaLambda0        = %"S_FLOAT"\n"

#define P_Omepnu    SET"OmegaNeutrino       = %"S_FLOAT"\n"
#define P_Omepsnu   SET"OmegaSterileNeutrino = %"S_FLOAT"\n"
#define P_Numnu     SET"Number of Neutrino  = %"S_FLOAT"\n"
#define P_NumMsvnu  SET"Number of Massive Neutrino  = %"S_INT"\n"
#define P_Hierachynu SET"Neutrino Hierarchy Type = %"S_INT"\n"

#define P_wlam0     SET"w0 of DE (CPL)      = %"S_FLOAT"\n"
#define P_wlam1     SET"wa of DE (CPL)      = %"S_FLOAT"\n"
#define P_Cs_DE     SET"Sound Speed of DE   = %"S_FLOAT"\n"

#define P_Omei      SET"OmegaMatterI        = %"S_FLOAT"\n"
#define P_fNL       SET"fNL                 = %"S_FLOAT"\n"
#define P_gNL       SET"gNL                 = %"S_FLOAT"\n"
#define P_Hubble    SET"Hubble              = %"S_FLOAT"\n"
#define P_nPS       SET"nPS                 = %"S_FLOAT"\n"
#define P_Boxsize   SET"Boxsize(Mpc/h)      = %"S_FLOAT"\n"

#define P_b8        SET"bias8               = %"S_FLOAT"\n"
#define P_bcmb      SET"biascmb             = %"S_FLOAT"\n"
#define P_Pkrescale SET"ampRescaling flag   = %"S_INT"\n"
#define P_As        SET"As(CMB,LCDM =2.076e-9) = %"S_FLOAT"\n"

#define P_Amax      SET"Amax                = %"S_FLOAT"\n"
#define P_Anow      SET"Anow                = %"S_FLOAT"\n"
#define P_Nx        SET"Nx                  = %"S_INT"\n"
#define P_Ny        SET"Ny                  = %"S_INT"\n"
#define P_Nz        SET"Nz                  = %"S_INT"\n"
#define P_Nspace    SET"Nspace              = %"S_INT"\n"
#define P_Mx        SET"Mx                  = %"S_LONG"\n"
#define P_My        SET"My                  = %"S_LONG"\n"
#define P_Mz        SET"Mz                  = %"S_LONG"\n"
#define P_MxMy      SET"MxMy                = %"S_LONG"\n"
#define P_Lx        SET"Lx                  = %"S_DOUBLE"\n"
#define P_Ly        SET"Ly                  = %"S_DOUBLE"\n"
#define P_Lz        SET"Lz                  = %"S_DOUBLE"\n"
#define P_Theta     SET"Theta               = %"S_FLOAT"\n"
#define P_Rsphere   SET"Rsphere             = %"S_FLOAT"\n"
#define P_Pradius   SET"Pradius             = %"S_FLOAT"\n"
#define P_ken       SET"Ken                 = %"S_FLOAT"\n"
#define P_ktot      SET"Ktot                = %"S_FLOAT"\n"
#define P_const0    SET"Const0              = %"S_FLOAT"\n"
#define P_poten0    SET"Poten0              = %"S_FLOAT"\n"
#define P_fact1     SET"Fact1               = %"S_FLOAT"\n"
#define P_fact2     SET"Fact2               = %"S_FLOAT"\n"
#define P_pfact     SET"Pfact               = %"S_FLOAT"\n"
#define P_stepnum   SET"Stepnum             = %"S_INT"\n"
#define P_nskip     SET"Nskip               = %"S_INT"\n"
#define P_local_nz  SET"Local_nz            = %"S_INT"\n"
#define P_local_sz  SET"Local_z_start       = %"S_INT"\n"
#define P_Zinit     SET"Z_init              = %"S_FLOAT"\n"
#define P_Rsmooth   SET"Rsmooth             = %"S_FLOAT"\n"
#define P_Rthooth   SET"Rth                 = %"S_FLOAT"\n"
#define P_Astep     SET"Astep               = %"S_FLOAT"\n"
#define P_Nstep     SET"Nstep               = %"S_INT"\n"
#define P_Istep     SET"Stepcount           = %"S_INT"\n"
#define P_Iseed     SET"Iseed               = %"S_INT"\n"
#define P_rvfile    SET"rvfilename          = %"S_STRING"\n"
#define P_rvprefix  SET"rvprefix            = %"S_STRING"\n"
#define P_Np        SET"Np                  = %"S_LONG"\n"
#define P_Zmin      SET"Zmin                = %"S_FLOAT"\n"
#define P_Zmax      SET"Zmax                = %"S_FLOAT"\n"
#define P_Powerflag SET"Powerflag           = %"S_INT"\n"
#define P_DEMODEL   SET"DarkEnergyModel     = %"S_STRING"\n"
#define P_Powfile   SET"Binary Powerfile    = %"S_STRING"\n"
#define P_Apowfile  SET"Ascii Powerfile     = %"S_STRING"\n"
#define P_XYZSHIFT  SET"XYZ_Shift_Flag      = %"S_INT"\n"
#define P_PTSIZE    SET"Particle_TYPE_Size  = %"S_INT"\n"
#define P_Viewer    SET"ViewerTrackFile     = %"S_STRING"\n"
#define P_Closing   "#End of Ascii Header\n"
#define P_NULL      "#------------------\n"
#define P_BLANK     "\n"

#define MPI_DATA(frw,wp,sp,simpar){\
 	ncnt += frw(wp,P_Myid,sp simpar.myid);\
	ncnt += frw(wp,P_Nid,sp simpar.nid);\
}
#define ENERGY(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_ken,sp simpar.ken);\
	ncnt += frw(wp,P_ktot,sp simpar.ktot);\
	ncnt += frw(wp,P_const0,sp simpar.const0);\
	ncnt += frw(wp,P_poten0,sp simpar.poten0);\
}
#if DEMODEL == VAR_DE
#define DE_MODEL(frw,wp,sp,simpar) {\
	ncnt += frw(wp,"### DE Model: PPF/FLUID/AXIONEFFECTIVEFLUID/EARLYQUINTESSENCE\n");\
	ncnt += frw(wp,"### As is the amplitude of P(k) at zi.\n");\
	ncnt += frw(wp,P_Omepnu,sp simpar.omepnu);\
	ncnt += frw(wp,P_Omepsnu,sp simpar.omepnusterile);\
    ncnt += frw(wp,P_DEMODEL   ,simpar.DarkEnergyModel);\
    ncnt += frw(wp,P_wlam0   ,sp simpar.wlam0);\
    ncnt += frw(wp,P_wlam1   ,sp simpar.wlam1);\
    ncnt += frw(wp,P_Cs_DE   ,sp simpar.CsDE);\
	ncnt += frw(wp,P_Pkrescale,sp simpar.Pk_rescaling);\
	if(simpar.Pk_rescaling ==0) ncnt += frw(wp,P_As,sp simpar.As);\
}
#else
#define DE_MODEL(frw,wp,sp,simpar) {\
    ncnt += frw(wp,P_wlam0   ,sp simpar.wlam0);\
    ncnt += frw(wp,P_wlam1   ,sp simpar.wlam1);\
    ncnt += frw(wp,P_Cs_DE   ,sp simpar.CsDE);\
	ncnt += frw(wp,"### bias8 is inverse of sigma_8 for Rsm=8 Mpc/h at z=0.\n");\
	ncnt += frw(wp,"### biascmb is inverse of sigma_8 for Rsm=130 Mpc/h at z=cmb.\n");\
	ncnt += frw(wp,"### Pk rescaling is the method to rescale power amp: \n");\
	ncnt += frw(wp,"###  -->  0(no rescaling)/1(@z=0)/2(@z=cmb)\n");\
	ncnt += frw(wp,"###  -->  When 1/2 is selected, bias8/biascmb should be provided.\n");\
	ncnt += frw(wp,"###  -->  If not, we set biascmb=1.\n");\
	ncnt += frw(wp,P_Pkrescale,sp simpar.Pk_rescaling);\
	if(simpar.Pk_rescaling ==1 ){\
		ncnt += frw(wp,P_b8,sp simpar.bias8);\
	}\
	else if(simpar.Pk_rescaling ==2 ){\
		ncnt += frw(wp,P_bcmb,sp simpar.biascmb);\
	}\
}
#endif

#define CORE(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### INITIAL CONDITION: 1 (ZA)/ 2 (2'nd order LPT)\n");\
 	ncnt += frw(wp,P_IC,sp simpar.IC);\
	ncnt += frw(wp,P_Omep,sp simpar.omep);\
	ncnt += frw(wp,P_Omepb,sp simpar.omepb);\
	ncnt += frw(wp,P_Omeplamb,sp simpar.omeplam);\
	ncnt += frw(wp,P_fNL,sp simpar.fNL);\
	ncnt += frw(wp,P_gNL,sp simpar.gNL);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### Hubble parameter is in 100km/sec/Mpc.\n");\
	ncnt += frw(wp,"### nPS is the power spectral index.\n");\
	ncnt += frw(wp,"### Boxsize should be in unit of Mpc/h.\n");\
	ncnt += frw(wp,P_Hubble,sp simpar.hubble);\
	ncnt += frw(wp,P_nPS,sp simpar.npow);\
	ncnt += frw(wp,P_Boxsize,sp simpar.boxsize);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	DE_MODEL(frw,wp,sp,simpar);\
	ncnt += frw(wp,"### powerflag = 0/1/2 (0: generate on the run, \n");\
	ncnt += frw(wp,"###  --> 1: read from camb file, 2: read from ascii file\n");\
	ncnt += frw(wp,P_Powerflag,sp simpar.powreadflag);\
	if(DEMODEL==1 || simpar.powreadflag==2) ncnt += frw(wp,P_Apowfile,   simpar.inpapkfilename);\
	else if(simpar.powreadflag==1){\
		ncnt += frw(wp,P_Powfile,   simpar.powfilename);\
	}\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Amax,sp simpar.amax);\
	ncnt += frw(wp,P_Anow,sp simpar.anow);\
	ncnt += frw(wp,P_Astep,sp simpar.astep);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,"### Nx, Ny, and Nz should be same & nspace should be int.\n");\
	ncnt += frw(wp,P_Nx,sp simpar.nx);\
	ncnt += frw(wp,P_Ny,sp simpar.ny);\
	ncnt += frw(wp,P_Nz,sp simpar.nz);\
	ncnt += frw(wp,P_Nspace,sp simpar.nspace);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Theta,sp simpar.theta);\
	ncnt += frw(wp,P_Nstep,sp simpar.nstep);\
	ncnt += frw(wp,P_Istep,sp simpar.stepcount);\
	ncnt += frw(wp,P_stepnum,sp simpar.stepnum);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Iseed,sp simpar.iseed);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_rvfile,   simpar.rvfilename);\
	ncnt += frw(wp,P_rvprefix,   simpar.rvprefix);\
	ncnt += frw(wp,P_BLANK);\
	ncnt += frw(wp,P_NULL);\
}
#define AUX1(frw,wp,sp,simpar){\
	ncnt += frw(wp,"### Miscellaneous parameters\n");\
	ncnt += frw(wp,P_Omei,sp simpar.omei);\
	ncnt += frw(wp,P_Mx,sp simpar.mx);\
	ncnt += frw(wp,P_My,sp simpar.my);\
	ncnt += frw(wp,P_Mz,sp simpar.mz);\
	ncnt += frw(wp,P_MxMy,sp simpar.mxmy);\
	ncnt += frw(wp,P_Lx,sp simpar.lnx);\
	ncnt += frw(wp,P_Ly,sp simpar.lny);\
	ncnt += frw(wp,P_Lz,sp simpar.lnz);\
}
#define DEFAULT_PARAMS(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_Starting);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_Closing);\
}

#define VIEWER(frw,wp,sp,simpar,ncnt){\
	ncnt += frw(wp,P_Viewer,    (simpar.Viewerfile));\
}
#define VIEWER2(frw,wp,sp,simpar,ncnt){\
	ncnt += frw(wp,P_Viewer,    (simpar->Viewerfile));\
}

#define PARAMS(frw,wp,sp,simpar){\
	ncnt += frw(wp,P_Starting);\
	CORE(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_XYZSHIFT,sp simpar.xyzshiftflag);\
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
	ncnt += frw(wp,P_XYZSHIFT,sp simpar.xyzshiftflag);\
	simpar.ptypesize = sizeof(pmparticletype);\
	ncnt += frw(wp,P_PTSIZE,sp simpar.ptypesize);\
	AUX1(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_Rsphere,sp simpar.sphere_radius);\
	ncnt += frw(wp,P_Pradius,sp simpar.particle_radius);\
	ENERGY(frw,wp,sp,simpar);\
	ncnt += frw(wp,P_fact1,sp simpar.fact1);\
	ncnt += frw(wp,P_fact2,sp simpar.fact2);\
	ncnt += frw(wp,P_pfact,sp simpar.pfact);\
	ncnt += frw(wp,P_nskip,sp simpar.nskip);\
	ncnt += frw(wp,P_Zmax,sp simpar.zmax);\
	ncnt += frw(wp,P_Zmin,sp simpar.zmin);\
	ncnt += frw(wp,P_local_nz,sp simpar.local_nz);\
	ncnt += frw(wp,P_local_sz,sp simpar.local_z_start);\
	ncnt += frw(wp,P_Np,sp simpar.np);\
	ncnt += frw(wp,P_NULL);\
	ncnt += frw(wp,P_Closing);\
}
