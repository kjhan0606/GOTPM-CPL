#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include <mpi.h>
#include <omp.h>
#include <unistd.h>
#include "pmheader.h"
#include "fof.h"
#include "Time.h"
#define MIN(a,b) a > b? b: a
#define MAX(a,b) a > b? a: b
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))
//particle p;
size_t readparticle(FoFTPtlStruct **,size_t ,int ,int ,int ,char *,int, MPI_Comm );
POSTYPE fof_link = 0.2;
POSTYPE lx,ly,lz;
float vscale,rscale;

#define TERMINAL 99

void xhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].x < ra[j+1].x) j++;
            if (rra.x < ra[j].x) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
void yhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].y < ra[j+1].y) j++;
            if (rra.y < ra[j].y) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}
void zhpsort(unsigned long n, particle ra[])
{
    unsigned long i,ir,j,l;
    particle rra;

    if (n < 2) return;
    l=(n >> 1)+1;
    ir=n;
    for (;;) {
        if (l > 1) {
            rra=ra[--l];
        } else {
            rra=ra[ir];
            ra[ir]=ra[1];
            if (--ir == 1) {
                ra[1]=rra;
                break;
            }
        }
        i=l;
        j=l+l;
        while (j <= ir) {
            if (j < ir && ra[j].z < ra[j+1].z) j++;
            if (rra.z < ra[j].z) {
                ra[i]=ra[j];
                i=j;
                j <<= 1;
            } else j=ir+1;
        }
        ra[i]=rra;
    }
}

/* (C) Copr. 1986-92 Numerical Recipes Software 71.+I0>+. */

HaloQ haloproperty(particle *member,size_t nmem){
	HaloQ halo;
	particle *pp;
	size_t i,j,k;
	double cx,cy,cz;
	double vx,vy,vz;
	POSTYPE xmin,xmax,ymin,ymax,zmin,zmax;
	POSTYPE x1,x2,y1,y2,z1,z2;
	halo.np=nmem;
	cx=cy=cz=vx=vy=vz=0;
	pp=member;
	xmax = ymax = zmax = -1.E25;
	xmin = ymin = zmin = +1.E25;
	for(i=0;i<nmem;i++){
		xmax = MAX(xmax,(pp->x));
		xmin = MIN(xmin,(pp->x));
		ymax = MAX(ymax,(pp->y));
		ymin = MIN(ymin,(pp->y));
		zmax = MAX(zmax,(pp->z));
		zmin = MIN(zmin,(pp->z));
		pp++;
	}
	if(xmin <= fof_link && xmax >=lx-fof_link){
		xhpsort(nmem,member-1);
		pp=member+1;
		for(i=1;i<nmem;i++){
			if(((pp->x)-(pp-1)->x) > 1.5*fof_link) {
				break;
			}
			pp++;
		}
		for(j=i;j<nmem;j++){
			member[j].x -= lx;
		}
	}
    if(ymin <= fof_link && ymax >=ly-fof_link){
                yhpsort(nmem,member-1);
                pp=member+1;
                for(i=1;i<nmem;i++){
                        if(((pp->y)-(pp-1)->y) > 1.5*fof_link) {
                                break;
                        }
						pp++;
                }
                for(j=i;j<nmem;j++){
                        member[j].y -= ly;
                }
     }
     if(zmin <= fof_link && zmax >=lz-fof_link){
                zhpsort(nmem,member-1);
                pp=member+1;
                for(i=1;i<nmem;i++){
                        if(((pp->z)-(pp-1)->z) > 1.5*fof_link) {
                                break;
                        }
						pp++;
                }
                for(j=i;j<nmem;j++){
                        member[j].z -= lz;
                }
     }


	pp=member;
	for(i=0;i<nmem;i++){
		cx+=(pp->x);
		cy+=(pp->y);
		cz+=(pp->z);
		vx+=pp->vx;
		vy+=pp->vy;
		vz+=pp->vz;
		pp++;
	}
	/* position in h^-1 Mpc */
	halo.x = cx/(POSTYPE)nmem *rscale;
	halo.y = cy/(POSTYPE)nmem *rscale;
	halo.z = cz/(POSTYPE)nmem *rscale;
	/* velocity in km/sec */
	halo.vx = vx/(float)nmem*vscale;
	halo.vy = vy/(float)nmem*vscale;
	halo.vz = vz/(float)nmem*vscale;
	return halo;
}
int myid,nid,mid;
int amyid, anid;
int main(int argc,char *argv[]){
	FILE *hfp,*pfp;
	double std,mean;
	int ntmp;
	size_t num;
	POSTYPE tmpx,tmpy,tmpz,dist2;
	POSTYPE fplmf,ptlmass;
	size_t i,j,k;
//	int si;
//	int nfof;
//	int iseed;
	int nowfile;
	int N,M;
	size_t N3;
	size_t np;
	size_t nhalo;
	int nstep;
	FoFTPtlStruct *ptl;
	FoFTPtlStruct *ptr;
	FoFBeginEndTree beginend;
	Box box;
	FoFTStruct *TREE;
	float wtime;
	int nfile;
	float size,hubble,npower,omep,omepb,omeplam,bias,smooth;
	int nx,ny,nz,nspace;
	float ntree,ntree1,theta;
	float zinit,amax,astep,anow,a;
	char halofile[190],memparticlefile[190];
	FILE *fp;
	POSTYPE zmin,zmax,zminlocal;
	size_t npadd,npwrite,npread;
	int initfile,finalfile;
	char infile[190];
	char hostname[128];
	long long ntreemax = 1000000L;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&amyid);
	MPI_Comm_size(MPI_COMM_WORLD,&anid);



	if(argc != 3) {
		if(myid==0){
			fprintf(stderr,"Error in # of arguments\n");
			fprintf(stderr,"%%fof fileleadingheadername nstep\n");
		}
		MPI_Finalize();
		exit(199);
	}

	MPI_Comm WorkWorld; 
	MPI_Group worldgroup,workgroup;

	{ 
		MPI_Comm_group(MPI_COMM_WORLD, &worldgroup);
		int nwork = anid-1; 
		int ranks[nwork];
		for(i=1;i<anid;i++) ranks[i-1] = i;

		MPI_Group_incl(worldgroup, nwork, ranks, &workgroup);

		if(amyid >0){
			MPI_Comm_create_group(MPI_COMM_WORLD, workgroup, 0,&WorkWorld);
			MPI_Comm_rank(WorkWorld, &myid);
			MPI_Comm_size(WorkWorld, &nid);
		}
	}
	printf("P%d has G%d\n", amyid, myid);




	FILE *ffp;
	double r2kineticfact,HSUB;
	SimParameters read_head(FILE *);
	nstep = atoi(argv[2]);
	sprintf(infile,"%s.%.5d%.5d",argv[1],nstep,amyid);
	if(amyid==0){
		ffp = fopen(infile,"r");
		simpar = read_head(ffp);
		fclose(ffp);
	}
	MPI_Bcast(&simpar,sizeof(SimParameters),MPI_BYTE,0,MPI_COMM_WORLD);
	nfile = simpar.nid;
	size = simpar.boxsize;
	hubble = simpar.hubble;
	npower = simpar.npow;
	omep = simpar.omep;
	omepb = simpar.omepb;
	omeplam = simpar.omeplam;
	bias = simpar.bias;
	smooth = simpar.rsmooth;
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	theta = simpar.theta;
	zinit = simpar.amax-1.;
	astep = simpar.astep;
	anow = simpar.anow;
       ny=nz=nx;
	a = anow;
	amax = simpar.amax;
	N3 = (size_t)(nx/nspace)*(size_t)(ny/nspace)*(size_t)(nz/nspace)/nfile;
	N3 = N3*0.5;
	rscale = size/nx;
	HSUB = sqrt(omep*pow3(amax/a)+omeplam+(1.-omep-omeplam)*pow2(amax/a));
	vscale = size/amax*a*a*100.*HSUB;

#ifndef MulOut
	sprintf(halofile,"FoF_halo_cat.%.5d",nstep);
	sprintf(memparticlefile,"FoF_member_particle.%.5d",nstep);
#else
	sprintf(halofile,"FoF_halo_cat.%.5d.%.5d",nstep, amyid);
	sprintf(memparticlefile,"FoF_member_particle.%.5d.%.5d",nstep, amyid);
#endif

	if(amyid==0) {
		hfp=fopen(halofile,"w");
		pfp=fopen(memparticlefile,"w");
		fwrite(&size,sizeof(float),1,hfp);
		fwrite(&hubble,sizeof(float),1,hfp);
		fwrite(&npower,sizeof(float),1,hfp);
		fwrite(&omep,sizeof(float),1,hfp);
		fwrite(&omepb,sizeof(float),1,hfp);
		fwrite(&omeplam,sizeof(float),1,hfp);
		fwrite(&bias,sizeof(float),1,hfp);
		fwrite(&nx,sizeof(int),1,hfp);
		fwrite(&nspace,sizeof(int),1,hfp);
		fwrite(&amax,sizeof(float),1,hfp);
		fwrite(&astep,sizeof(float),1,hfp);
		fwrite(&anow,sizeof(float),1,hfp);
		fclose(hfp);
		fclose(pfp);

		printf("P%d size = %g hubble = %g\n",amyid,size,hubble);
		printf("P%d npower = %g omep = %g omeplam = %g bias = %g smooth=%g\n",
				amyid,npower,omep,omeplam,bias,smooth);
		printf("P%d zinit = %g astep = %g anow = %g \n",amyid,zinit,astep,anow);
		printf("P%d nx = %d ny= %d nz= %d nspace=%d\n",amyid,nx,ny,nz,nspace);
		printf("P%d rscale = %g vscale= %g\n",amyid,rscale,vscale);
		int icount = anid-1;
		int isyncRead = 0;
		do{
			int src,ok,tag;
            MPI_Status rstatus;

	        MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG,MPI_COMM_WORLD,&rstatus);

			while(isyncRead >= maxReads && rstatus.MPI_TAG == 2){
	        	MPI_Probe(MPI_ANY_SOURCE, 3,MPI_COMM_WORLD,&rstatus);
			}
            src = rstatus.MPI_SOURCE;
			tag = rstatus.MPI_TAG;
			if(tag==3){
            	MPI_Recv(&ok,1,MPI_INT,src,tag,MPI_COMM_WORLD,&rstatus);
				isyncRead --;
				printf("P%d <- G%d: a report of the end of reading (isync= %d/%d)\n", 
						myid, src-1, isyncRead, maxReads);
			}
			else if(tag==2){
            	MPI_Recv(&ok,1,MPI_INT,src,tag,MPI_COMM_WORLD,&rstatus);
            	MPI_Send(&ok,1,MPI_INT,src,tag,MPI_COMM_WORLD);
				isyncRead ++;
				printf("P%d -> G%d: a permission to read (isync=%d/%d)\n", 
						myid, src-1, isyncRead, maxReads);
			}
			else if(tag ==1){
            	MPI_Recv(&ok, 1, MPI_INT,src, tag, MPI_COMM_WORLD,&rstatus);
				MPI_Send(&ok, 1, MPI_INT,src, 2, MPI_COMM_WORLD);
            	MPI_Recv(&ok, 1, MPI_INT,src, 3, MPI_COMM_WORLD,&rstatus);
				printf("P%d <- G%d: a request to dump halo data\n", myid, src-1);
			}
			else if(tag == TERMINAL){
            	MPI_Recv(&ok, 1, MPI_INT,src, tag, MPI_COMM_WORLD,&rstatus);
				icount --;
				printf("P0 <- G%d: a report to finalize the whole job (icountdown=%d)\n", src-1, icount);
			}
			fflush(stdout);
		}while(icount>0);

	}
	else {
#ifdef MulOut
		hfp=fopen(halofile,"w");
		pfp=fopen(memparticlefile,"w");
		fclose(hfp);
		fclose(pfp);
#endif
          
		fof_link = fof_link * nspace;


		gethostname(hostname,127);
		if((ptl = (FoFTPtlStruct *) malloc(sizeof(FoFTPtlStruct)*100)) == NULL){
			fprintf(stderr,"Error allocating ptl %ld\n",N3);
			exit(99);
		}
		M = 1;
		if((TREE = (FoFTStruct *) malloc(sizeof(FoFTStruct)*ntreemax)) == NULL){ 
			fprintf(stderr,"Error allocating TREE\n"); 
			exit(99); 
		}
	
		wtime = WALLCLOCK();
		(void)WALLCLOCK();
		box.x = box.y = box.z = 0.;
		box.width = nx;
		lx = nx;
		ly = ny;
		lz = nz;
		printf("G%d has Lx=%g Ly=%g Lz=%g\n",myid,lx,ly,lz);
		np = 0;
		npwrite = 0;
		/* realintifile & realfinalfile are the real file range
		 * but (initfile, finalfile) could be a fake-file range to be set to synchronize the 
		 * parallel process */
		int realinitfile,realfinalfile;
		int LocalNFiles = (nfile-1 +nid)/nid;
		int n_pes = (nfile + LocalNFiles - 1)/LocalNFiles;
		{
			if(myid >= n_pes) {
				/* These are for a fake  to redundant ranks to whom no files are allocated */
				initfile = nfile;
				finalfile = nfile + LocalNFiles;
				realinitfile = realfinalfile = -1;
			}
			else {
				realinitfile = myid * LocalNFiles;
				if(myid == n_pes-1) realfinalfile = nfile;
				else realfinalfile = realinitfile + LocalNFiles;
				initfile = realinitfile;
				finalfile = initfile + LocalNFiles;
			}
		}
		int nthreads;
#pragma omp parallel
		{
#pragma omp master
			{
				nthreads = omp_get_num_threads();
			}
		}
		particle *linked = (particle *) malloc(sizeof(particle)*MaxLinkedParticles*nthreads);
		INDXTYPE *linkedtindx = (INDXTYPE *) malloc(sizeof(INDXTYPE)*MaxLinkedParticles*nthreads);
		for(nowfile=initfile;nowfile<finalfile;nowfile++){
			sprintf(infile,"%s.%.5d%.5d",argv[1],nstep,nowfile);
			npadd = readparticle(&ptl,np,nstep,nowfile,nz,infile,N3, WorkWorld);
			np += npadd;
	
			if(nowfile == finalfile-1){
//				MPI_Barrier(WorkWorld);
				int src,dest;
				src = (myid+1+nid)%nid;
				dest = (myid-1+nid)%nid;
				if(myid == n_pes-1) { src = 0;}
				if(myid == 0) { dest = n_pes-1;}
				if(myid < n_pes) {
					if(1){
						char file[190];
						sprintf(file,"BottomFaceContactHalo%.5d%.5d.dat",src,nstep);
						FILE *ffp = fopen(file,"r"); fseek(ffp, 0L, SEEK_END);
						long fsize = ftell(ffp);
						npread = fsize/sizeof(particle);
						fclose(ffp);
					}
					else{
						MPI_Status status;
						MPI_Sendrecv(&npwrite,sizeof(size_t),MPI_BYTE,dest,0,
								&npread,sizeof(size_t),MPI_BYTE,src,0,WorkWorld,&status);
					}
					ptl = (FoFTPtlStruct*)realloc(ptl,
							sizeof(FoFTPtlStruct)*(np+npread));
					ReadBottomFaceContact(ptl+np,npread,linked,src,nstep,nz);
					np += npread;
				}
			}
			else {
				ptl = (FoFTPtlStruct *)realloc(ptl,sizeof(FoFTPtlStruct)*np);
			}
			HaloBound *halobound;
			if((nowfile < nfile || nowfile == finalfile-1) && np >0)
			{
				zmin = 2.E23;
				zmax = -2.E23;
#pragma omp parallel for reduction(min:zmin) reduction(max:zmax)
				for(i=0;i<np;i++){
					zmin = MIN(zmin,ptl[i].r[2]);
					zmax = MAX(zmax,ptl[i].r[2]);
				}
		
				if(nowfile==initfile) zminlocal = zmin;
	
	
				ptl = (FoFTPtlStruct *)realloc(ptl,sizeof(FoFTPtlStruct)*np);
				if(np /MIN_CELL_PARTICLE_NUM*2.0  > ntreemax){
					ntreemax = np/MIN_CELL_PARTICLE_NUM*2.0;
					TREE = (FoFTStruct *)realloc(TREE,sizeof(FoFTStruct)*ntreemax);
				}
				int nbin = 2*nthreads;
				double nsize = (double)nx/(double)nbin;
				double gridSize = (double)size/(double)nbin;
				printf("G%d has zmin/zmax = %g %g & omp pixel size = %g for nthreads= %d with ntreemax= %ld\n",
						myid, zmin, zmax, gridSize, nthreads, ntreemax);

#pragma omp parallel for
				for(i=0;i<np;i++){
					ptl[i].type = TYPE_PTL;
					ptl[i].sibling = &ptl[i+1];
					ptl[i].gridLL = NULL;
					ptl[i].haloindx = -1;
					ptl[i].tindx = i;
				}
				ptl[np-1].sibling = NULL;
#ifndef _OPENMP
				FoF_Make_Tree(TREE,ptl,np,box);
#else
				void Omp_FoF_Make_Tree(FoFTStruct *,FoFTPtlStruct *,size_t,Box,int,long long );
				Omp_FoF_Make_Tree(TREE,ptl,np,box, nx,ntreemax);
#endif
				printf("G%d has made fof tree for nbin= %d\n", myid,nbin);fflush(stdout);
	
				if(gridSize < 10.){
					fprintf(stderr,"Error: small pixelsize %g h^-1 Mpc nthreads= %d nx= %d\n", gridSize, nthreads, nx);
					MPI_Finalize();
					return 1;
				}
				FoFTPtlStruct **Grids1D = (FoFTPtlStruct **)malloc(sizeof(FoFTPtlStruct*)*nbin);
				size_t *npGrids1D = (size_t *)malloc(sizeof(size_t)*nbin);
				for(i=0;i<nbin;i++) {
					Grids1D[i] = NULL;
					npGrids1D[i] = 0;
				}
#pragma omp parallel private(i)
				{
					int idthread = omp_get_thread_num();
					for(i=0;i<np;i++){
						int ix = ptl[i].r[0]/nsize;
						if(ix >= nbin) ix = nbin-1;
						if( (ix/2) == idthread){
							FoFTPtlStruct *tp = Grids1D[ix];
							Grids1D[ix] =  ptl + i;
							ptl[i].gridLL = tp;
							npGrids1D[ix] ++;
						}
					}
				}
//				size_t halostartid[nbin], halofinalid[nbin];
				size_t *halostartid = (size_t*)malloc(sizeof(size_t)*nbin);
				size_t *halofinalid = (size_t*)malloc(sizeof(size_t)*nbin);
	
				printf("G%d passed the linked list\n", myid);
				/* First half */
				for(i=0;i<2;i++) {
					/*
					if(i==0) printf("G%d(%s) Before first half\n", myid, hostname);
					else printf("G%d(%s) Before second half\n", myid, hostname);
					*/
#pragma omp parallel 
					{
						POSTYPE llx,lly,llz;
						llx = lx; lly = ly; llz = lz;
						int idthread = omp_get_thread_num();
						particle *threadLinked = linked+MaxLinkedParticles*idthread;
						INDXTYPE *threadLinkedTindx = linkedtindx+MaxLinkedParticles*idthread;
						int ix = idthread*2 + i;
						size_t haloid=0;
						size_t ii;
						for(ii=0;ii<ix;ii++) haloid += npGrids1D[ii];
						// this is the temporary starting point of halo id 
						halostartid[ix] = haloid; 
						FoFTPtlStruct *target = Grids1D[ix];
						while(target){
							if(target->included == NO){
								particle pp;
								pp.x = target->r[0];
								pp.y = target->r[1];
								pp.z = target->r[2];
								size_t nmem=pnew_fof_link(&pp,fof_link,TREE,ptl,
										threadLinked,threadLinkedTindx,
										haloid,llx,lly,llz);
		
								if(nmem >= MinNumMem){
									haloid ++;
								}
								else{
									POSTYPE hzmin,hzmax;
									hzmin = 2*nz;
									hzmax = -nz;
									for(ii=0;ii<nmem;ii++){
										hzmin = MIN(hzmin, threadLinked[ii].z);
										hzmax = MAX(hzmax, threadLinked[ii].z);
									}
									if(hzmin < zminlocal+fof_link || hzmax > zmax-fof_link) haloid ++;
									else {
										for(ii=0;ii<nmem;ii++){
											ptl[threadLinkedTindx[ii]].haloindx = -1;
										}
									}
								}
							}
							target = target->gridLL;
						}
						// this is temporary end point of halo id 
						halofinalid[ix] = haloid;
//						printf("G%d with T%d finishes\n", myid, idthread);
					}
					if(i==0) printf("G%d(%s) passed the first half\n", myid, hostname);
					else if(i==1) printf("G%d(%s) passed the second half\n", myid, hostname);
				}
				printf("G%d(%s) End of halo finding\n", myid, hostname);fflush(stdout);
				long long *halostack = (long long *)malloc(sizeof(long long)*(halofinalid[nbin-1]));
				/* stacking the halo id to the left (lower value) */
				/* halostack has the map from the current haloid to a new stacked haloid */
				size_t istack = 0;
				for(i=0;i<nbin;i++){
					size_t istart = halostartid[i];
					size_t ifinal = halofinalid[i];
					size_t jj;
					for(jj=istart;jj<ifinal;jj++){
						halostack[jj] = istack;
						istack ++;
					}
				}

				if(0){
					nhalo = 0;
					for(i=0;i<np;i++){
						long long jj = ptl[i].haloindx;
						if(jj >= 0) {
							ptl[i].haloindx = halostack[jj]; 
							nhalo = MAX(nhalo, halostack[jj]);
						}
					}
					/* increase the number of halo by one for the c convention*/
					nhalo ++;

				}
				else{
					nhalo = istack;
					long idmax,idmin;
					idmax = -10000;
					idmin = 100000000L;
					for(i=0;i<np;i++){
						long long jj = ptl[i].haloindx;
						if(jj>=0){
							ptl[i].haloindx = halostack[jj];
							idmax = MAX(idmax,halostack[jj]);
							idmin = MIN(idmin,halostack[jj]);
						}
					}
					printf("G%d has halo idmin/max = %ld %ld for nhalo=%ld\n", myid,idmin,idmax,nhalo);
					if(nhalo <= idmax){
						printf("ERRORORROROROR in halo id %ld %ld\n", idmax,nhalo);
					}
				}


				free(halostack);
				free(halofinalid); free(halostartid); 
				free(npGrids1D); free(Grids1D);
	
				/*
				for(i=0;i<np;i++){
					if(ptl[i].included == NO){
						p.x = ptl[i].r[0];
						p.y = ptl[i].r[1];
						p.z = ptl[i].r[2];
						num=pnew_fof_link(&p,fof_link,TREE,ptl,linked,nhalo,lx,ly,lz);
						nhalo ++;
					}
				}
				*/
				printf("G%d has %ld halos from %ld particles\n",myid,nhalo,np);
				halobound = NULL;
				halobound = (HaloBound *)malloc(sizeof(HaloBound)*nhalo);
				if(halobound == NULL){
					printf("error in allocating halobound for %ld\n", nhalo);
				}
				CheckHaloBound(nhalo,halobound,ptl,np,fof_link,
					zmin,zmax, zminlocal);
				printf("G%d has checked halo boundary\n",myid);fflush(stdout);
			}
			{
				if(nowfile != finalfile-1) {
					int flag;
					WriteIsolatedHalo(nhalo,halobound,ptl,linked,halofile,
							memparticlefile);
					printf("G%d has saved isolated halos\n",myid);fflush(stdout);
					if(nowfile == initfile) flag = 0;
					else flag = 1;
					npwrite += WriteBottomFaceContact(nhalo,halobound,ptl,linked,
							flag,nstep, myid);
#ifndef PREVIOUS_VERSION
					if(nowfile < nfile)
#endif
					{
						np = StackUpContactParticleLeftWard(nhalo,halobound,ptl,np);
						printf("G%d has saved Bottom Faced Halo %ld\n",myid,npwrite);fflush(stdout);
						free(halobound);
					}
				}
				if(nowfile == finalfile-1){
					WriteFinalHalo(nhalo,halobound,ptl,linked,halofile,
							memparticlefile);
					int terminal = -1;
					int rtag = TERMINAL;
					MPI_Send(&terminal, 1, MPI_INT,0,rtag, MPI_COMM_WORLD);
					MPI_Finalize();
					return 0;
				}
			}
		}
	} // End of amyid !=0
	if(amyid==0){
		printf("############################################\n");
		printf("Congratulations!!!!\n");
		printf("\n\n Now the end of the this opfof.exe\n\n");
		printf("############################################\n");
	}
	MPI_Finalize();
}
