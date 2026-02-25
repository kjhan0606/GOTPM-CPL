#include<stdio.h>
#include<stddef.h>
#include<stdlib.h>
#include<string.h>
#include <mpi.h>
#include <math.h>
#include "pmheader.h"
#include "params.h"
void determine_mpi_misc_param(SimParameters *simpar){
	int nid,myid;
	float zi;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	simpar->nid = nid;
	simpar->myid = myid;
#ifdef XYZDBL
	simpar->xyzshiftflag = 1;
#else
	simpar->xyzshiftflag = 0;
#endif
	zi = simpar->amax-1;
	simpar->omei = simpar->omep*pow(1+zi,3)/(simpar->omep*pow(1+zi,3) +
			simpar->omeplam*pow(1+zi,3*(1+simpar->wlam0))+(1-simpar->omep-simpar->omeplam)*pow(1+zi,2));
	simpar->mx = simpar->nx/simpar->nspace;
	simpar->my = simpar->ny/simpar->nspace;
	simpar->mz = simpar->nz/simpar->nspace;
	simpar->mxmy = simpar->mx * simpar->my;
	simpar->lnx = simpar->nx;
	simpar->lny = simpar->ny;
	simpar->lnz = simpar->nz;
	simpar->sphere_radius = 4;
	simpar->particle_radius = 4;
	simpar->rth = 8/(simpar->boxsize/simpar->nx);
	simpar->zinit = simpar->amax -1;
}
