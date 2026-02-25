#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>


#include "pmheader.h"
#include "pm_common.h"
#include "flow.h"
#include "varpm.h"
#include "Memory.h"
#include "fortran.h"
#include "fftw3.h"
#include "fftw3-mpi.h"

#if DEMODEL == 0
void poissoncork(float *fftwmesh){
	long i, npix;
	npix = (long)simpar.local_nz*(long)simpar.ny*(long)simpar.nx;
	for(i=0;i<npix;i++) fftwmesh[i] *= simpar.pcorr;
}

#elif DEMODEL == 1
void poissoncork(float *fftwmesh){
	long i,j,k;
	long nxh = (simpar.nx/2+1);
	long mx = 2*(simpar.nx/2+1);
	double wj,wk,pcork;
	float wave;
	for(j=0;j<simpar.local_nz;j++){
		wj = (j + simpar.local_z_start);
		for(k=0;k<simpar.nz;k++){
			wk = k;
			for(i=0;i<=nxh;i++){
				wave = sqrt(i*i + wj*wj + wk*wk);
				pcork = interppcork(wave);
				fftwmesh[2*i   + mx*(k + simpar.nz*j)] *= pcork;
				fftwmesh[2*i+1 + mx*(k + simpar.nz*j)] *= pcork;
			}
		}
	}
}
#endif



void pmmain(int *nowflag, int *now){
	float *fftwmesh,poten,efold;
	int nzheight,nxpos,nsizepmparticletype;
	int iflagwholeden,iflagbinnedden;
	float pmas;
#ifdef INCLUDE_TREE_FORCE
    efold=0.9;
#else
    efold = 0.5;
#endif
    nxpos = (int)offsetof(pmparticletype,x) /sizeof(float);
    nsizepmparticletype = (int)sizeof(pmparticletype) /sizeof(float);
	{
		long tlnp,lnp;
		double rngc;
		rngc = (float)simpar.nx*(float)simpar.ny*(float)simpar.nz;
		lnp = (long ) np;
		MPI_Reduce(&lnp,&tlnp,1,MPI_LONG,MPI_SUM,0,MPI_COMM_WORLD);
		MPI_Bcast(&tlnp,1,MPI_LONG,0,MPI_COMM_WORLD);
		pmas = rngc/(double)tlnp;

	}

	{
        nzheight = ceil(simpar.zmax) - (int)simpar.zmin;
        printf("P%d has np=%ld %d from %g to %g and local_z_start=%d simpar.local_nz=%d\n",
                simpar.myid,np,nzheight,simpar.zmin,simpar.zmax,simpar.local_z_start,simpar.local_nz);
        den = (float *) Malloc(((nzheight+4)*(2*(simpar.nx/2+1))*simpar.ny)*
                sizeof(float), PPTR(den));
        if(simpar.local_nz > 0) fftwmesh = (float *) Malloc(((simpar.local_nz+0)*(2*(simpar.nx/2+1))*simpar.ny)*
                sizeof(float),PPTR(fftwmesh));
        else fftwmesh = (float *) Malloc(sizeof(float),PPTR(fftwmesh));

        VarPM_tsc2(den,simpar.nx,simpar.ny,simpar.nz, simpar.zmin,simpar.zmax, pmparticles,
                nsizepmparticletype,nxpos, np,pmas,
                fftwmesh, simpar.local_nz,simpar.local_z_start);
        {/* Saving x-z slice file */
            Savexzslice(fftwmesh);
        }

        if(halfstep.first != 1 || halfstep.second !=1 ) {
            if((iflagwholeden=flagwholeden(simpar.amax,simpar.anow,simpar.astep))){
                SaveWholeDen(fftwmesh);
            }
            if((iflagbinnedden=flagBinnedDen(simpar.amax,simpar.anow,simpar.astep))){
                SaveBinnedData(fftwmesh);
            }
        }
        /* This is forward FFT */
//        fftforward_(fftwmesh);
		TIMER_START(11);
		void fftwndtransforward_(float *);
        fftwndtransforward_(fftwmesh);
		TIMER_STOP(11);
		if(simpar.myid==0) fprintf(stdout,"CPU(F-FFTW3)   = %f\n",ELAPSED_TIME(11));
		/* */
//		poissoncork(fftwmesh);
// pcorr test - minseong 2024_06_15

        /* This is to measure the correlation function */
        if(flagpsmeasure(simpar.amax,simpar.anow,simpar.astep,simpar.stepcount)){
            void correl_(float *,float *,int *,int *,float *,int *,
                    float *,float *,float *,float *,float *,float *);
            float *pk; int nstep;
            pk = (float *)Malloc(sizeof(float)*6*(simpar.nx+1),PPTR(pk));
            nstep = simpar.nx+1;
            correl_(fftwmesh,&simpar.anow,nowflag,now,&simpar.rth,&simpar.nstep,
                    pk,pk+nstep,pk+2*nstep,pk+3*nstep,pk+4*nstep,pk+5*nstep);
            Free(pk);
        }
        psolver_(fftwmesh,&poten,&efold);
//        fftbackward_(fftwmesh);
		TIMER_START(11);
        void fftwndtransbackward_(fftw_complex *);
        fftwndtransbackward_((fftw_complex*) fftwmesh);
		TIMER_STOP(11);
		if(simpar.myid==0) fprintf(stdout,"CPU(B-FFTW3)   = %f\n",ELAPSED_TIME(11));
		if(0){
			void GottGrid(float *);
			if(simpar.stepcount %10 == 1)GottGrid(fftwmesh);
        /* This is for the four point FDA */
		}
        {
            int nzfloor,izslice,izslicep1,local_nz_varpm;
            simpar.zmin = simpar.zmin; simpar.zmax = simpar.zmax;
            /* No buffer slab is needed */
            den = (float*)Realloc(den,nzheight*(2*(simpar.nx/2+1))*simpar.ny*sizeof(float));
            adVarPM_fda4_mesh(fftwmesh,simpar.local_nz,simpar.local_z_start,den,simpar.zmin,
                    simpar.zmax, simpar.nx,simpar.ny,simpar.nz);
            Free(fftwmesh);
            nzfloor = (int) simpar.zmin;
            local_nz_varpm = nzheight;

            den = (float*)Realloc(den,sizeof(float)*2*(simpar.nx/2+1)*simpar.ny*(nzheight+7));
            /* This is for the padding of boundary slices */
            paddingpotentialmesh(den,nzfloor,nzheight,local_nz_varpm,simpar.nx,simpar.ny,simpar.nz,-3, 4);
            /* WARNING local_nz --> local_nz_varpm*/
            {
                void fda4(pmparticletype *,int,float *,int,int,int,int);
                fda4(pmparticles,np,den,local_nz_varpm,(int)simpar.zmin,halfstep.first,halfstep.second);
            }
        }
        Free(den);
	}
}
