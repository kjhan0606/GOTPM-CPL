#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
#include"Memory.h"
#include"tree.h"
#include"force_spline.h"
#define READY 0
#define WRITING 1
#define NP_TAG 111
#define R_TAG 456
#define MAX(a,b) ((a)>(b)?(a):(b))
#define MIN(a,b) ((a)<(b)?(a):(b))
#define pow2(a) ((a)*(a))
#define pow3(a) ((a)*(a)*(a))

FILE *outwp;

typedef long indxtype;

typedef struct HaloInfo{
	float mass,cx,cy,cz;
	int indx,mpeak,isubhalo;
	float spin,sigv,rotL,radius,q,s;
	float v[3][3],ang[3],pvx,pvy,pvz;
	float Es; /* This is the surface pressure term defined in Shaw 2006 */
}HaloInfo;
#define MAXNH 60000000
HaloInfo *ahalo,shaloinfo,rhaloinfo;
int nhalo;
typedef struct vector3d{
	float x,y,z;
}vector3d;
vector3d *r,*vr;
indxtype  *indx;
float d[3],v[3][3];
double cx,cy,cz;
double cvx,cvy,cvz;
double ang_mom[3];
double tenergy,sim2real,para_spin,rotation,tmass;
double mass,pvr,ampang;
double sigv;
int motherrank=0;
float radius;
float q,s,r1[3],r2[3],r3[3];

float real2sim,massscale,rscale;


int sorthalo(void *a, void *b){
	HaloInfo *aa;
	HaloInfo *bb;
	aa = (HaloInfo*)a;
	bb = (HaloInfo*)b;
	if(aa->indx > bb->indx)
		return 1;
	else if(aa->indx < bb->indx)
		return -1;
	else if(aa->isubhalo >bb->isubhalo)
		return 1;
	else if(aa->isubhalo <bb->isubhalo)
		return -1;
	else return 0;
}

int myid,nid;
int shalonum,rhalonum;
int mpeak;
int ready=READY,writing=WRITING;
int snp,rnp;
/*                              */
float amax,a,rng,size,hubble;
int ng,nspace;
float omep,omepb,omeplam,epsilon;
int nx,ny,nz;
/*                              */
double onesolarmass=1.989E33L;
double com2real,real2com,potentfact;
double pntmass,r1kineticfact,r2kineticfact;
double G = 6.672E-8L;
double pc = 3.08567802E18L;
double pi = 3.14159265354979L;
float redshift,rng,omepk,Hsub;


void Save_data(HaloInfo);
void triaxialshape_(int *,float *,float*,float*,float *,float *,float *,float *,double *,double*,double*);

int physical_parameters(void);
void pspline_(void);
void i_potent_spline(void);

/*
#include "force_spline.h"
*/
void halo_potent(int nend,vector3d *r,float *penergy){
    TStruct *TREE;
    TPtlStruct *ptl;
    particle p;
    Box box;
    float theta=1.;
    int i,j,k;
    static int np;
    np = nend;
    TREE = (TStruct *) Malloc(sizeof(TStruct)*np*4,PPTR(TREE));
    ptl = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*np,PPTR(ptl));
    for(i=0;i<np;i++){
        ptl[i].type = TYPE_PTL;
        ptl[i].r[0] = r[i].x; ptl[i].r[1] = r[i].y; ptl[i].r[2] = r[i].z;
        ptl[i].mass = 1;
    }

    box = findbox(ptl,np);
#ifdef DEBUG
		fprintf(stderr,"Now before make_tree for np= %d box x,y,z = %g %g %g width = %g\n",np,
				box.x,box.y,box.z,box.width);
		if(np==23){
			for(i=0;i<np;i++){
				fprintf(stderr,"%g %g %g : %g %g %g\n",r[i].x,r[i].y,r[i].z,ptl[i].r[0],ptl[i].r[1],ptl[i].r[2]);
			}
		}
#endif
    Make_Tree(TREE,ptl,np,box);
#ifdef DEBUG
		fprintf(stderr,"Now after make_tree \n");
#endif
    for(i=0;i<np;i++){
        p.x = ptl[i].r[0]; p.y = ptl[i].r[1]; p.z = ptl[i].r[2];
        penergy[i] = treeplumpotential(&p,theta,TREE,ptl)*potentfact;
        
    }
    Free(ptl);
    Free(TREE);
}
void halo_properties(int nend,vector3d *r,vector3d *vr){
    float distx,disty,distz,tmpvx,tmpvy,tmpvz;
    float distvx,distvy,distvz;
    float tmpvvx,tmpvvy,tmpvvz;
    float *kenergy,*penergy,*rv;
    float absdist,rvx,rvy,rvz;
    float *x,*y,*z;
    double mrv;
	float potmin;
    static int np;
    int i,j,k;
    np = nend;

    penergy = (float *) Malloc(sizeof(float)*np,PPTR(penergy));
    halo_potent(np,r,penergy);
#ifdef MAXDEN
	{
		int Numnear = 27;
		Numnear = MIN(Numnear,np);
		DenPeakCenter(r,np,&cx,&cy,&cz,Numnear);
	}
#elif MAXPOT
	{
		potmin = 0.;
		for(i=0;i<np;i++){
			if(penergy[i]<potmin){
				j = i;
				potmin = penergy[i];
			}
		}
		cx = r[j].x; cy = r[j].y; cz = r[j].z;
	}
#elif COM
	{
		cx = cy = cz=0;
		for(i=0;i<np;i++){
			cx += r[i].x; cy += r[i].y; cz += r[i].z;
		}
		cx = cx/(double)np; cy = cy/(double)np; cz = cz/(double)np;
	}
#else
	{
		fprintf(stderr,"Error.. Please choose one of MAXDEN, MAXPOT, and  COM \n");
	}

#endif



    cvx = cvy = cvz = 0.;
    for(i=0;i<np;i++){
        cvx += vr[i].x; cvy += vr[i].y; cvz += vr[i].z;
    }
    cvx = cvx/ (float) np; cvy = cvy/(float) np; cvz = cvz/(float) np;
    kenergy = (float *) Malloc(sizeof(float)*np,PPTR(kenergy));
    rv = (float *) Malloc(sizeof(float)*np,PPTR(rv));
	sigv = 0.L;
    for(i=0;i<np;i++){
        distx = (r[i].x-cx)*r1kineticfact;
        disty = (r[i].y-cy)*r1kineticfact;
        distz = (r[i].z-cz)*r1kineticfact;
        tmpvx = (vr[i].x-cvx)*r2kineticfact;
        tmpvy = (vr[i].y-cvy)*r2kineticfact;
        tmpvz = (vr[i].z-cvz)*r2kineticfact;
        tmpvvx = (distx+tmpvx); tmpvvy = (disty+tmpvy); tmpvvz = (distz+tmpvz);
        kenergy[i] = 0.5*(tmpvvx*tmpvvx+tmpvvy*tmpvvy+tmpvvz*tmpvvz);
		sigv += kenergy[i]+kenergy[i];
    }
	sigv = sqrt(sigv/(double)np)/1.E5L;

    for(i=0;i<3;i++){
        ang_mom[i] = 0.L;
    } 
    mrv = 0.L;
    for(i=0;i<np;i++){
        distx = r[i].x - cx; disty = r[i].y - cy; distz = r[i].z - cz;
        distvx = vr[i].x-cvx; distvy = vr[i].y-cvy; distvz = vr[i].z-cvz;
        ang_mom[0] += disty*distvz-distz*distvy;
        ang_mom[1] += distz*distvx-distx*distvz;
        ang_mom[2] += distx*distvy-disty*distvx;
        absdist = sqrt(distx*distx+disty*disty+distz*distz);
        rvx = disty*distvz-distz*distvy;
        rvy = distz*distvx-distx*distvz;
        rvz = distx*distvy-disty*distvx;
        rv[i] = sqrt(rvx*rvx+rvy*rvy+rvz*rvz)/absdist;
        mrv = MAX(mrv,rv[i]);
    }
    tenergy = 0.;
    for(i=0;i<np;i++) tenergy += penergy[i]*0.5 + kenergy[i];
    para_spin = sqrt(sqr(ang_mom[0])+sqr(ang_mom[1])+sqr(ang_mom[2])) *
        sqrt(-tenergy)/G/pntmass*sim2real*
        r2kineticfact/pow((double)np,2.5);
    /* angular momentum in the unit of M_odot kpc km/sec */
    rotation = sqrt(sqr(ang_mom[0])+sqr(ang_mom[1])+sqr(ang_mom[2])) *
        sim2real/(1.E3L*pc)*r2kineticfact*pntmass/onesolarmass;
    tmass = pntmass*np;
    {
        x = (float *) Malloc(sizeof(float)*np,PPTR(x));
        y = (float *) Malloc(sizeof(float)*np,PPTR(y));
        z = (float *) Malloc(sizeof(float)*np,PPTR(z));
        for(i=0;i<np;i++){
            x[i] = r[i].x; y[i] = r[i].y; z[i] = r[i].z;
        }
#ifdef DEBUG
		fprintf(stderr,"Now before triaxialshape cx,y,z = %g %g %g\n",cx,cy,cz);
#endif
        triaxialshape_(&np,x,y,z,&q,&s,&d[0],&v[0][0],&cx,&cy,&cz);
#ifdef DEBUG
		fprintf(stderr,"Now after triaxialshape\n");
#endif
        Free(x);Free(y);Free(z);
    }
    Free(rv);Free(penergy);Free(kenergy);
}
void getparameter(FILE *wp){
    float bias,rsmooth;
    float npow0;
    float z_init,delta_a;
    float SPHERE_RADIUS,particle_radius,rtheta;
    fscanf(wp,"%f %f\n",&size,&hubble);
    fscanf(wp,"%f %f %f %f %f %f\n",&npow0,&omep,&omepb,&omeplam,&bias,&rsmooth);
    fscanf(wp,"%d %d %d %d\n",&nx,&ny,&nz, &nspace);
    fscanf(wp,"%f %f %f\n",&SPHERE_RADIUS,&particle_radius,&rtheta);
    fscanf(wp,"%f %f %f\n",&z_init,&delta_a,&a);
    rng = nx;
    amax = 1.+z_init;
    redshift = amax/a - 1.;
}
typedef struct bparticle{
	float x,y,z,vx,vy,vz;
	float dist;
} bparticle;
int sortdist(const void *a,const void *b){ 
	double aa,bb; 
	aa = ((bparticle*)a)->dist; 
	bb = ((bparticle*)b)->dist;
	if(aa < bb) return -1; 
	else if(aa > bb) return +1; 
	else return 0;
}
void getEs(int np, vector3d *r, vector3d *vr, HaloInfo *halo){
	int i,j,k;
	int i80,i90;
	double tmpx,tmpy,tmpz,res;
	double distx,disty,distz;
	double tmpvx,tmpvy,tmpvz;
	double tmpvvx,tmpvvy,tmpvvz;
	bparticle *bp;

	i80 = np*0.8;
	i90 = np*0.9;
	if(i80>=np-1 || i90>=np-1) {
		halo->Es = 0.;
	}
	else {
		bp = (bparticle *)Malloc(sizeof(bparticle)*np,PPTR(bp));
		for(i=0;i<np;i++){
			bp[i].x = r[i].x;
			bp[i].y = r[i].y;
			bp[i].z = r[i].z;
			bp[i].vx = vr[i].x;
			bp[i].vy = vr[i].y;
			bp[i].vz = vr[i].z;
			tmpx = bp[i].x - cx;
			tmpy = bp[i].y - cy;
			tmpz = bp[i].z - cz;
			bp[i].dist = sqrt(tmpx*tmpx+tmpy*tmpy+tmpz*tmpz);
		}
		qsort(bp,np,sizeof(bparticle),sortdist);
		res = 0.L;
		for(i=i80;i<np;i++){
   	    	distx = (bp[i].x-cx)*r1kineticfact;
   	     	disty = (bp[i].y-cy)*r1kineticfact;
	   	    distz = (bp[i].z-cz)*r1kineticfact;
   		    tmpvx = (bp[i].vx-cvx)*r2kineticfact;
   		    tmpvy = (bp[i].vy-cvy)*r2kineticfact;
   		    tmpvz = (bp[i].vz-cvz)*r2kineticfact;
   		    tmpvvx = (distx+tmpvx); tmpvvy = (disty+tmpvy); tmpvvz = (distz+tmpvz);
   		    res += (tmpvvx*tmpvvx+tmpvvy*tmpvvy+tmpvvz*tmpvvz);/*Kinetic energy per unit particle mass */
		}
		halo->Es = pow(bp[i90].dist,3.L)/(pow(bp[np-1].dist,3.L)-pow(bp[i80].dist,3.L))*res;
		Free(bp);
	}
}
float find_r(vector3d*,int,float,float,float*,float*,float*,double,double,double);
void GethaloInfo(vector3d *r,vector3d *vr,int np, HaloInfo *halo){
	int i,j,k;
	float rfact,vfact;
	rfact = rng/size;
	vfact = 1.E5/r2kineticfact;
	for(i=0;i<np;i++){
		r[i].x = r[i].x *rfact;
		r[i].y = r[i].y *rfact;
		r[i].z = r[i].z *rfact;
		vr[i].x = vr[i].x *vfact;
		vr[i].y = vr[i].y *vfact;
		vr[i].z = vr[i].z *vfact;
	}
	{
		/*
		int Numnear = 27;
		Numnear = MIN(Numnear,np);
		DenPeakCenter(r,np,&cx,&cy,&cz,Numnear);
		*/
	}
#ifdef DEBUG
	fprintf(stderr,"Now before halo properties\n");
#endif
	{
		halo_properties(np,r,vr);
	}
#ifdef DEBUG
	fprintf(stderr,"Now after halo properties\n");
#endif
	{
        float *px,*py,*pz;
        px = (float *) Malloc(sizeof(float)*10000,PPTR(px));
        py = (float *) Malloc(sizeof(float)*10000,PPTR(py));
        pz = (float *) Malloc(sizeof(float)*10000,PPTR(pz));
        r1[0] = v[0][0]; r1[1] = v[0][1]; r1[2] = v[0][2];
        r2[0] = v[1][0]; r2[1] = v[1][1]; r2[2] = v[1][2];
        r3[0] = v[2][0]; r3[1] = v[2][1]; r3[2] = v[2][2];
        radius = find_r(r,np,q,s,r1,r2,r3,cx,cy,cz)*size/rng;
        Free(pz);
        Free(py);
        Free(px);
	}
#ifdef DEBUG
	fprintf(stderr,"Now after find_r\n");
#endif

    {
        cvx = cvy = cvz = 0;
        for(j=0;j<np;j++){ 
			cvx += vr[j].x; 
			cvy += vr[j].y; 
			cvz += vr[j].z; 
		}
        cvx = cvx/(double)np;
        cvy = cvy/(double)np;
        cvz = cvz/(double)np;
		/*
        sigv = 0.L;
        for(j=0;j<np;j++){
			double distx,disty,distz;
			double tmpx,tmpy,tmpz,tx,ty,tz;
			distx = r1kineticfact*(r[j].x-cx); 
			disty = r1kineticfact*(r[j].y-cy); 
			distz = r1kineticfact*(r[j].z-cz);
			tmpx = r2kineticfact*(vr[j].x-cvx);
			tmpy = r2kineticfact*(vr[j].y-cvy);
			tmpz = r2kineticfact*(vr[j].z-cvz);
			tx = distx+tmpx; ty = disty+tmpy; tz = distz+tmpz;
			sigv += tx*tx+ty*ty+tz*tz;
        }
        sigv = sigv/(double) np;
        sigv = sqrt(sigv);
        sigv = sigv/1.E5L;
		*/
    }
    mass = np*pntmass/onesolarmass;
    pvr = sqrt(cvx*cvx+cvy*cvy+cvz*cvz)*r2kineticfact/1.E5;
    if(isnan(para_spin) != 0) {
    	para_spin = -9999.;
    }
    {
        ampang = pow2(ang_mom[0])+pow2(ang_mom[1])+pow2(ang_mom[2]);
        ampang = sqrt(ampang);
        ang_mom[0] = ang_mom[0]/ampang;
        ang_mom[1] = ang_mom[1]/ampang;
        ang_mom[2] = ang_mom[2]/ampang;
    }
#ifdef DEBUG
	fprintf(stderr,"Now before getEs\n");
#endif
	{
		halo->cx = cx*size/rng;
		halo->cy = cy*size/rng;
		halo->cz = cz*size/rng;
		halo->q = q;
		halo->s = s;
		halo->spin = para_spin;
		halo->pvx = cvx*r2kineticfact/1.E5;
		halo->pvy = cvy*r2kineticfact/1.E5;
		halo->pvz = cvz*r2kineticfact/1.E5;
		halo->radius = radius;
		halo->sigv = sigv;
		halo->rotL = rotation;
		halo->ang[0] = ang_mom[0];
		halo->ang[1] = ang_mom[1];
		halo->ang[2] = ang_mom[2];
		halo->mass = mass;
		for(i=0;i<3;i++)for(j=0;j<3;j++)halo->v[i][j] = v[i][j];
		/*
		getEs(np, r,vr,halo);
		*/
	}
#ifdef DEBUG
	fprintf(stderr,"Now after getEs\n");
#endif
}



/*
int main(int argc, char *argv[]) {
*/
int MAIN_(int argc, char *argv[]){
	MPI_Status mstatus,cstatus;
	MPI_Request request;
	char filetag[80],filename[80],parafile[190];
	int i,j,snp,rnp,flag;
	int finish=-99999;
	int snend=100,rnend;
	int nstep,index;
	int src,dest;
	float tmp_tidal,tmp;
	FILE *fp,*paraid;
	char rheader[80],rfilename[80],wfilename[80];



	(void )MPI_Init(&argc,&argv);
	(void )MPI_Comm_size(MPI_COMM_WORLD,&nid);
	(void )MPI_Comm_rank(MPI_COMM_WORLD,&myid);

	if(argc != 2){
		if(myid == motherrank)
		fprintf(stderr,".... mpimf.exe [nstep]\n");
		exit(1);
	}
	nstep = atoi(argv[1]);

	if(Make_Total_Memory()==0){
		fprintf(stderr,"P%d Error initializing mem %ldMB\n",myid,NMEG);
	}

	i_potent_spline();
	if(myid==motherrank){
		ahalo = (HaloInfo*)Malloc(sizeof(HaloInfo)*MAXNH,PPTR(ahalo));
		sprintf(parafile,"params.%.5d",nstep);
		paraid = fopen(parafile,"r");
		getparameter(paraid);fclose(paraid);
	}
	MPI_Bcast(&ng,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&nx,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&ny,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&nz,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&nspace,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omep,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&rng,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omepb,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omeplam,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&size,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&hubble,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&amax,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&a,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&redshift,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);

	omepk = 1.-omep-omeplam;
   	Hsub = sqrt(omep*pow3(amax/a)+omeplam+omepk*pow2(amax/a));
   	r1kineticfact = size/hubble/rng/amax*a*100.E5*hubble*Hsub;
   	r2kineticfact = size/hubble/amax*a*a*100.E5*hubble*Hsub;
   	sim2real = size/hubble/rng/amax*a*1.E6*pc;
   	real2sim = 1./sim2real;
   	pntmass = 3./8./pi/G*sqr(100.E5*hubble)*1.E6*pc;
   	pntmass = pntmass*pow3(size/hubble)/pow3(rng/nspace)*omep;
   	potentfact = G*pntmass/sim2real;
   	massscale = pntmass/onesolarmass*hubble;
   	rscale = size/nx;
	/*
   	if(myid==motherrank) printf("masscale=%g rscale = %g\n",massscale,rscale);
	*/
	printf("P%d has omep=%g rng=%g r1kineticfact=%g r2kineticfact=%g potentfact=%g\n",
			myid,omep,rng,r1kineticfact,r2kineticfact,potentfact);
	i_potent_spline();
	if(myid==motherrank){
		char outfile[190];
		sprintf(outfile,"n%.5dhalomass.dat",nstep);
		outwp=fopen(outfile,"w");
	}

	sprintf(rfilename,"HFIND.DATA.%.5d",nstep);
	if(myid==motherrank) {
		fprintf(stderr,"Opening %s file for input\n",rfilename);
	}
	shalonum = rhalonum = 0;
	if(myid == motherrank) {
		index = 0;
		fp = fopen(rfilename,"r");
		while(fread(&mpeak,sizeof(int),1,fp) == 1){
			for(i=0;i<mpeak;i++){
				fread(&snp,sizeof(int),1,fp);
				r = (vector3d *) Malloc(sizeof(vector3d)*snp,PPTR(r));
				vr = (vector3d *) Malloc(sizeof(vector3d)*snp,PPTR(vr));
				indx = (indxtype*)Malloc(sizeof(indxtype)*snp,PPTR(indx));
				fread(r,sizeof(vector3d),snp,fp);
				fread(vr,sizeof(vector3d),snp,fp);
				fread(indx,sizeof(indxtype),snp,fp);
				if(snp <5) {
					Free(indx);Free(vr);Free(r);
					continue;
				}
				do {
					MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
					src = mstatus.MPI_SOURCE;
					dest = mstatus.MPI_SOURCE;
					MPI_Recv(&ready,1,MPI_INT,src,READY,MPI_COMM_WORLD,&cstatus);
					if(ready == READY){
						shalonum++;
						MPI_Send(&snp,1,MPI_INT,dest,NP_TAG, MPI_COMM_WORLD);
						MPI_Send(r,snp*sizeof(vector3d),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
						MPI_Send(vr,snp*sizeof(vector3d),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
						MPI_Send(indx,snp*sizeof(indxtype),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
						MPI_Send(&index,1,MPI_INT,dest,R_TAG,MPI_COMM_WORLD);
						MPI_Send(&mpeak,1,MPI_INT,dest,R_TAG,MPI_COMM_WORLD);
						MPI_Send(&i,1,MPI_INT,dest,R_TAG,MPI_COMM_WORLD);
#ifdef DEBUG
						printf("P%d sending halo data to P%d with np=%d for %d\n",myid,dest,snp,index);
#endif
					}
					else {
						rhalonum++;
						MPI_Probe(src,NP_TAG,MPI_COMM_WORLD,&cstatus);
						MPI_Recv(&rhaloinfo,sizeof(HaloInfo),MPI_BYTE,src,NP_TAG,MPI_COMM_WORLD,&cstatus);
#ifdef DEBUG
						printf("P%d receiving halo data from P%d\n",myid,src);
#endif

						Save_data(rhaloinfo);
					}
				} while(ready != READY); 
				Free(indx);Free(vr);Free(r);
			}
			index ++;
		}
		fclose(fp);
exithere :
		printf("complete................... \n");
		fflush(stdout);
		j = 0;
		for(i=rhalonum+1;i<=shalonum;){
			MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
			MPI_Recv(&ready,1,MPI_INT,mstatus.MPI_SOURCE,READY,
					MPI_COMM_WORLD,&cstatus);
			if(ready == WRITING){
				i++;
				MPI_Probe(mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD,&cstatus);
				MPI_Recv(&rhaloinfo,sizeof(HaloInfo),MPI_BYTE,mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD, &cstatus);
				Save_data(rhaloinfo);
			}
			else {
				j++;
				MPI_Send(&finish,1,MPI_INT,mstatus.MPI_SOURCE, NP_TAG,MPI_COMM_WORLD);
			}
		}
		for(i=1;i<nid-j;i++){
			MPI_Recv(&ready,1,MPI_INT,MPI_ANY_SOURCE,READY, MPI_COMM_WORLD,&mstatus);
			MPI_Send(&finish,1,MPI_INT,mstatus.MPI_SOURCE, NP_TAG,MPI_COMM_WORLD);
		}
	}
	/* 다른 프로세서는 여기서 계산을 한다  */
	else {
		ready = READY;
		MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
		MPI_Recv(&snp,1,MPI_INT,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
		while(snp != finish){
			long numstack;
			snend = MAX(snend,snp);
			r = Malloc(sizeof(vector3d)*snend,PPTR(r));
			vr = Malloc(sizeof(vector3d)*snend,PPTR(vr));
			indx = Malloc(sizeof(indx)*snend,PPTR(indx));
			MPI_Recv(r,sizeof(vector3d)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(vr,sizeof(vector3d)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(indx,sizeof(indxtype)*snp,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(&index,1,MPI_INT,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(&mpeak,1,MPI_INT,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(&i,1,MPI_INT,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
#ifdef DEBUG
			fprintf(stderr,"P%d received data from P%d with np=%d \n",myid,motherrank,snp);
#endif

			shaloinfo.indx = index;
			shaloinfo.mpeak = mpeak;
			shaloinfo.isubhalo = i;
			numstack = CurMemStack();
			GethaloInfo(r,vr,snp,&shaloinfo);
			InitialOldMemStack(numstack);
#ifdef DEBUG
			fprintf(stderr,"P%d sending data to P%d with np=%d \n",myid,motherrank,snp);
#endif

			ready = WRITING;
			MPI_Issend(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD, &request);
			MPI_Wait(&request,&mstatus);
			MPI_Send(&shaloinfo,sizeof(HaloInfo),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD);

			ready = READY;
			MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
			MPI_Recv(&snp,1,MPI_INT,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);

			Free(indx);Free(vr);Free(r);
		}
	}
	if(myid == motherrank) {
		HaloInfo *tt;
		fprintf(stderr,"End of calculation. Writing data\n");
		qsort(ahalo,nhalo,sizeof(HaloInfo),sorthalo);
		tt = ahalo;
		for(i=0;i<nhalo;i++){
			fprintf(outwp,"%g %f %f %f %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
					tt->mass,tt->cx,tt->cy,tt->cz,tt->indx,tt->mpeak,tt->isubhalo,tt->spin,tt->sigv,
					tt->rotL,tt->radius,tt->q,tt->s,tt->v[0][0],tt->v[0][1],tt->v[0][2],
					tt->v[1][0],tt->v[1][1],tt->v[1][2],tt->v[2][0],tt->v[2][1],tt->v[2][2],
					tt->ang[0],tt->ang[1],tt->ang[2],tt->pvx,tt->pvy,tt->pvz, tt->Es);
			tt++;

		}
		fclose(outwp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;
}
#define WRITEDOWNNOW {\
	int i,ii,inow,ic;\
	HaloInfo *tt;\
	qsort(ahalo,nhalo,sizeof(HaloInfo),sorthalo);\
	for(ii=0;ii<nhalo;ii++){\
		if(ahalo[ii+1].indx > ahalo[ii].indx+1){\
			break;\
		}\
	}\
	tt = ahalo;\
	for(i=0;i<ii;i++){\
		fprintf(outwp,"%g %f %f %f %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",\
				tt->mass,tt->cx,tt->cy,tt->cz,tt->indx,tt->mpeak,tt->isubhalo,tt->spin,tt->sigv,\
				tt->rotL,tt->radius,tt->q,tt->s,tt->v[0][0],tt->v[0][1],tt->v[0][2],\
				tt->v[1][0],tt->v[1][1],tt->v[1][2],tt->v[2][0],tt->v[2][1],tt->v[2][2],\
				tt->ang[0],tt->ang[1],tt->ang[2],tt->pvx,tt->pvy,tt->pvz, tt->Es);\
		tt++;\
	}\
	ic = 0;\
	for(i=ii;i<nhalo;i++){\
		ahalo[ic++] = ahalo[i];\
	}\
	nhalo = nhalo - ii;\
}
void Save_data(HaloInfo shalo){
	ahalo[nhalo++] = shalo;
	if(nhalo > 10000){
		WRITEDOWNNOW;
	}
}
