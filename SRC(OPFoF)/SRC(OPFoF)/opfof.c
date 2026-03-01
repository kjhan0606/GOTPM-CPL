#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include "mpi.h"
#include "fof.h"
#include "Time.h"
#include "pmheader.h"
#define MIN(a,b) a > b? b: a
#define MAX(a,b) a > b? a: b
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))
particle p;
size_t readparticle(FoFTPtlStruct **,size_t ,char *, char);
float fof_link = 0.2;
float lx,ly,lz;

	float size,hubble,npower,omep,omepb,omepl,bias,smooth;
	int nx,ny,nz,nspace;
	float ntree,ntree1,theta;
	float zinit,amax,astep,anow,a;
	int iseed;
	float vscale,pscale;
	FILE *hfp,*pfp;
	char halofile[190],memparticlefile[190];
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
	float xmin,xmax,ymin,ymax,zmin,zmax;
	float x1,x2,y1,y2,z1,z2;
	halo.np=nmem;
	cx=cy=cz=vx=vy=vz=0;
	pp=member;
	xmax = ymax = zmax = -1.E25;
	xmin = ymin = zmin = +1.E25;
	for(i=0;i<nmem;i++){
		xmax = MAX(xmax,pp->x);
		xmin = MIN(xmin,pp->x);
		ymax = MAX(ymax,pp->y);
		ymin = MIN(ymin,pp->y);
		zmax = MAX(zmax,pp->z);
		zmin = MIN(zmin,pp->z);
		pp++;
	}
	if(xmin <= fof_link && xmax >=lx-fof_link){
		xhpsort(nmem,member-1);
		pp=member+1;
		for(i=1;i<nmem;i++){
			if((pp->x-(pp-1)->x) > 1.5*fof_link) {
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
                        if((pp->y-(pp-1)->y) > 1.5*fof_link) {
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
                        if((pp->z-(pp-1)->z) > 1.5*fof_link) {
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
		cx+=pp->x;
		cy+=pp->y;
		cz+=pp->z;
#ifndef ONLYPOS
		vx+=pp->vx;
		vy+=pp->vy;
		vz+=pp->vz;
#endif
		pp++;
	}
	/* position in h^-1 Mpc */
	halo.x = cx/(float)nmem *pscale;
	halo.y = cy/(float)nmem *pscale;
	halo.z = cz/(float)nmem *pscale;
	/* velocity in km/sec */
#ifndef ONLYPOS
	halo.vx = vx/(float)nmem*vscale;
	halo.vy = vy/(float)nmem*vscale;
	halo.vz = vz/(float)nmem*vscale;
#endif
	return halo;
}
int myid,nid;
int main(int argc,char *argv[]){
	double std,mean;
	int ntmp;
	size_t num,ii;
	float tmpx,tmpy,tmpz,dist2;
	float fplmf,ptlmass;
	size_t i,j,k;
	int si,snp;
	int nowfile;
	int N,M;
	size_t N3;
	size_t np,addhere;
	int nstep;
	FoFTPtlStruct *ptl;
	FoFTPtlStruct *ptr;
	particle *linked;
	int nfof;
	FoFBeginEndTree beginend;
	Box box;
	FoFTStruct *TREE;
	float wtime;
	int nfile;
	FILE *fp;
	size_t nhalo;
	float zmin,zmax,zminlocal;
	size_t npadd,npwrite,npread;
	MPI_Status status;
	int initfile,finalfile;
	char pflag;

	MPI_Init(&argc,&argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);



	if(argc != 4) {
		fprintf(stderr,"Error in # of arguments\n");
		fprintf(stderr,"%%fof nstep nfile [O/N: old or new style of parameter]\n");
		exit(199);
	}
	else {
		FILE *ffp;
		char parafile[190];
		double r2kineticfact,HSUB;
		nstep = atoi(argv[1]);
		nfile = atoi(argv[2]);
		sprintf(parafile,"params.%.5d",nstep);
		if((ffp=fopen(parafile,"r"))==NULL){
			fprintf(stderr,"Error opening %s\n",parafile);
			exit(99);
		};
		sprintf(&pflag,"%s",argv[3]);
		if(pflag == 'O')  {
			fscanf(ffp,"%g %g",&size,&hubble);
			fscanf(ffp,"%g %g %g %g %g %g",&npower,&omep,&omepb,&omepl,&bias,&smooth);
			fscanf(ffp,"%d %d %d %d",&nx,&ny,&nz,&nspace);
			fscanf(ffp,"%g %g %g",&ntree,&ntree1,&theta);
			fscanf(ffp,"%g %g %g",&zinit,&astep,&anow);
		}
		else if(pflag == 'N'){
			SimParameters  read_sim_parameter_file(FILE*),simpar;
			simpar = read_sim_parameter_file(ffp);
			size = simpar.boxsize;
			hubble = simpar.hubble;
			npower = simpar.npow;
			omep = simpar.omep;
			omepb = simpar.omepb;
			omepl = simpar.omeplam;
			bias = simpar.bias;
			nx = simpar.nx;
			ny = simpar.ny;
			nz = simpar.nz;
			nspace = simpar.nspace;
			anow = simpar.anow;
			zinit = simpar.zinit;
			amax = simpar.amax;
		}
		else {
			fprintf(stderr,"Error in the N/O flag\n");
			exit(999);
		}
        ny=nz=nx;
		a = anow;
		fclose(ffp);
		N3 = (size_t)(nx/nspace)*(size_t)(ny/nspace)*(size_t)(nz/nspace)/nid;
		N3 = N3*1.5;
		pscale = size/nx;
		HSUB = sqrt(omep*pow3(amax/a)+omepl+(1.-omep-omepl)*pow2(amax/a));
		vscale = size/amax*a*a*100.*HSUB;
           
		fof_link = fof_link * nspace;
		printf("%g %g\n",size,hubble);
		printf("%g %g %g %g %g\n",npower,omep,omepl,bias,smooth);
		printf("%d %d %d %d\n",nx,ny,nz,nspace);
	}
	if((ptl = (FoFTPtlStruct *) malloc(sizeof(FoFTPtlStruct)*100)) == NULL){
		fprintf(stderr,"Error allocating ptl %ld\n",N3);
		exit(99);
	}
	M = 1;
	/*
	if((p = (particle *) malloc(sizeof(particle)*M)) == NULL){
                fprintf(stderr,"Error allocating p\n");
                exit(99);
        }
		*/
	linked = (particle *) malloc(sizeof(particle)*MaxLinkedParticles);
	if((TREE = (FoFTStruct *) malloc(sizeof(FoFTStruct)*
					N3/NODE_HAVE_PARTICLE*3)) == NULL){
                fprintf(stderr,"Error allocating TREE\n");
                exit(99);
        }
	wtime = WALLCLOCK();
	(void)WALLCLOCK();
	box.x = box.y = box.z = 0.;
	box.width = nx;
	sprintf(halofile,"FoF_halo_cat.%.5d",nstep);
	sprintf(memparticlefile,"FoF_member_particle.%.5d",nstep);
	if(myid==0){
		hfp=fopen(halofile,"w");
		pfp=fopen(memparticlefile,"w");
		fwrite(&size,sizeof(float),1,hfp);
		fwrite(&hubble,sizeof(float),1,hfp);
		fwrite(&npower,sizeof(float),1,hfp);
		fwrite(&omep,sizeof(float),1,hfp);
		fwrite(&omepb,sizeof(float),1,hfp);
		fwrite(&omepl,sizeof(float),1,hfp);
		fwrite(&bias,sizeof(float),1,hfp);
		fwrite(&nx,sizeof(int),1,hfp);
		fwrite(&nspace,sizeof(int),1,hfp);
		fwrite(&amax,sizeof(float),1,hfp);
		fwrite(&astep,sizeof(float),1,hfp);
		fwrite(&anow,sizeof(float),1,hfp);
		fclose(hfp);
	}
	lx = nx;
	ly = ny;
	lz = nz;
	printf("P%d has Lx=%g Ly=%g Lz=%g\n",myid,lx,ly,lz);
	np = 0;
	npwrite = 0;
	initfile = nfile/nid*myid;
	finalfile = nfile/nid*(myid+1);
	printf("P%d has a file range %d %d\n",myid,initfile,finalfile);
	for(nowfile=initfile;nowfile<finalfile;nowfile++){
		char infile[190];
		sprintf(infile,"SyncINITIAL.%.5d%.5d",nstep,nowfile);
		printf("P%d now opening %s\n",myid,infile);
		npadd = readparticle(&ptl,np,infile,pflag);
		np += npadd;
		if(nowfile == finalfile-1){
			int src,dest;
			src = (myid+1+nid)%nid;
			dest = (myid-1+nid)%nid;
			MPI_Sendrecv(&npwrite,sizeof(size_t),MPI_BYTE,dest,0,
					&npread,sizeof(size_t),MPI_BYTE,src,0,MPI_COMM_WORLD,&status);
			{
				ptl = (FoFTPtlStruct*)realloc(ptl,
						sizeof(FoFTPtlStruct)*(np+npread));
				ReadBottomFaceContact(ptl+np,npread,linked,src,nstep);
				np += npread;
			}
		}
		zmin = 2.E23;
		zmax = -2.E23;
		for(i=0;i<np;i++){
			zmin = MIN(zmin,ptl[i].r[2]);
			zmax = MAX(zmax,ptl[i].r[2]);
		}

		if(nowfile==initfile) zminlocal = zmin;

		printf("Now we have zmin=%g zmax=%g\n",zmin,zmax);
		ptl = (FoFTPtlStruct *)realloc(ptl,sizeof(FoFTPtlStruct)*np);
		for(i=0;i<np;i++){
			ptl[i].type = TYPE_PTL;
			ptl[i].sibling = &ptl[i+1];
			ptl[i].haloindx = -1;
		}
		ptl[np-1].sibling = NULL;
		FoF_Make_Tree(TREE,ptl,np,box);
		/* Local FoF and tag with a new halo number nhalo */
		nhalo = 0;
		for(i=0;i<np;i++){
			if(ptl[i].included == NO){
				p.x = ptl[i].r[0];
				p.y = ptl[i].r[1];
				p.z = ptl[i].r[2];
				num=pnew_fof_link(&p,fof_link,TREE,ptl,linked,nhalo,lx,ly,lz);
				nhalo ++;
			}
		}
		printf("P%d has %ld halos from %ld particles\n",myid,nhalo,np);
		{
			HaloBound *halobound;
			int flag;
			halobound = (HaloBound *)malloc(sizeof(HaloBound)*nhalo);
			CheckHaloBound(nhalo,halobound,ptl,np,fof_link,
				zmin,zmax, zminlocal);
			printf("P%d has checked halo boundary\n",myid);fflush(stdout);
			if(nowfile != finalfile-1) {
				WriteIsolatedHalo(nhalo,halobound,ptl,linked,halofile,
						memparticlefile);
			}
			else if(nowfile == finalfile-1){
				WriteFinalHalo(nhalo,halobound,ptl,linked,halofile,
						memparticlefile);
				return 0;
			}
			printf("P%d has saved isolated halos\n",myid);fflush(stdout);
			if(nowfile == initfile) flag = 0;
			else flag = 1;
			npwrite += WriteBottomFaceContact(nhalo,halobound,ptl,linked,
					flag,nstep);
			printf("P%d has saved Bottom Faced Halo %ld\n",myid,npwrite);fflush(stdout);
			np = StackUpContactParticleLeftWard(nhalo,halobound,ptl,np);
			free(halobound);
		}
	}
}
