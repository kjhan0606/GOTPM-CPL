/* 
 * #ifdef OLD is the memory devouring job.
 * */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include<mpi.h>





#include "Memory.h"
#include "merger.h"
#define READY 0
#define WRITING 1
#define NP_TAG 111
#define R_TAG 456
int ready=READY,writing=WRITING;
int myid,nid;




typedef struct PtlPos{
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
} PtlPos;
typedef struct PtlVel{
	float vx,vy,vz;
} PtlVel;






double eratio,ke,pe;

#define MIN(A,B) ((A)<(B) ? (A): (B))
#define MAX(A,B) ((A)>(B) ? (A): (B))
#define sqr(A) ((A)*(A))
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define sqr(A) ((A)*(A))
double sim2real,real2sim;
double pntmass,r1kineticfact,r2kineticfact;
double potentfact;
double onesolarmass=1.989E33L;
double  G= 6.672E-8L;
double pc= 3.08567802E18L;
double pi= 3.1415926535L;
double rotation;
int nspace;
float r1[3],r2[3],r3[3],radius;
float a;
float amax;
float redshift;
float size;
float rng,npower,bias;
float omep,omepb;
float omeplam, wlam0, wlam1;
float hubble;
float omepk,Hsub;
float astep,anow,pscale;
int nx,ny,nz;
double vscale;
float *den;
void findden(int , PtlPos *, int ,int ,int );
void triaxialshape_(int *,float *,float *,float *,float *,
		float *,float *,float *,double *,double *,double *);
float wc(float );
float find_r(PtlPos *, int, float ,float, float *, float *,float *);
double ang_mom[3];
float q,s;
float d[3],v[3][3];
double para_spin;
float tenergy,tmass,pvx;
double cx,cy,cz,cvx,cvy,cvz,pvr;
double sigv;
size_t Fread(void *buff,size_t size,size_t nmem,FILE *fp){
	char a,b,c,d;
	char *data;
	int i,j;
	size_t result;
	result = fread(buff,size,nmem,fp);
	data = (char *) buff;
	for(i=0;i<nmem*size;i+=4){
		a = data[i]; b = data[i+1]; c = data[i+2]; d = data[i+3];
		data[i] = d;
		data[i+1] = c;
		data[i+2] = b;
		data[i+3] = a;
	}
	return result;
}
/*
#define fread(a,b,c,d) Fread(a,b,c,d)
*/
#include "tree.h"
#include "Time.h"
#include "force_spline.h"
void halo_potent(int nend,PtlPos *r,float *penergy){
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
	/*
	{
		float xtmpmax,ytmpmax,ztmpmax;
		box.x = box.y =box.z =box.width =  -1.E20;
		for(i=0;i<np;i++){
			box.x = MIN(box.x,r[i].x);
			box.y = MIN(box.y,r[i].y);
			box.z = MIN(box.z,r[i].z);
			xtmpmax = MAX(xtmpmax,r[i].x);
			ytmpmax = MAX(ytmpmax,r[i].y);
			ztmpmax = MAX(ztmpmax,r[i].z);
		}
		box.width = MAX(box.width,xtmpmax-box.x);
		box.width = MAX(box.width,ytmpmax-box.y);
		box.width = MAX(box.width,ztmpmax-box.z);
	}
	box.x = box.y = box.z = 0.;box.width = rng;
	*/
	Make_Tree(TREE,ptl,np,box);
	for(i=0;i<np;i++){
		p.x = ptl[i].r[0]; p.y = ptl[i].r[1]; p.z = ptl[i].r[2];
		penergy[i] = treeplumpotential(&p,theta,TREE,ptl)*potentfact;
	}
	Free(ptl);
	Free(TREE);
}
INDXTYPE  halo_properties_kw(int nend,PtlPos *r,PtlVel *vr,INDXTYPE *indx,double *KW){
	float distx,disty,distz,tmpvx,tmpvy,tmpvz;
	float distvx,distvy,distvz;
	float tmpvvx,tmpvvy,tmpvvz;
	float *kenergy,*penergy,*rv;
	float absdist,rvx,rvy,rvz;
	float *x,*y,*z;
	double mrv;
	static int np;
	int i,j,k;
	INDXTYPE mbp;
	np = nend;
	cx = cy = cz = cvx = cvy = cvz = 0.;
	for(i=0;i<np;i++){
		cx += r[i].x; cy += r[i].y; cz += r[i].z;
		cvx += vr[i].vx; cvy += vr[i].vy; cvz += vr[i].vz;
	}
	cx = cx/(float) np; cy = cy/(float) np; cz = cz/(float) np;
	cvx = cvx/ (float) np; cvy = cvy/(float) np; cvz = cvz/(float) np;
	kenergy = (float *) Malloc(sizeof(float)*np,PPTR(kenergy));
	penergy = (float *) Malloc(sizeof(float)*np,PPTR(penergy));
	rv = (float *) Malloc(sizeof(float)*np,PPTR(rv));
	for(i=0;i<np;i++){
		distx = (r[i].x-cx)*r1kineticfact;
		disty = (r[i].y-cy)*r1kineticfact;
		distz = (r[i].z-cz)*r1kineticfact;
		tmpvx = (vr[i].vx-cvx)*r2kineticfact;
		tmpvy = (vr[i].vy-cvy)*r2kineticfact;
		tmpvz = (vr[i].vz-cvz)*r2kineticfact;
		tmpvvx = (distx+tmpvx); tmpvvy = (disty+tmpvy); tmpvvz = (distz+tmpvz);
		kenergy[i] = 0.5*(tmpvvx*tmpvvx+tmpvvy*tmpvvy+tmpvvz*tmpvvz);
	}
	{
		float xmin,ymin,zmin,xmax,ymax,zmax;
		xmin = ymin = zmin = 1.E20;
		xmax = ymax = zmax = -1.E20;
		for(i=0;i<np;i++){
			xmin = MIN(xmin,r[i].x);
			ymin = MIN(ymin,r[i].y);
			zmin = MIN(zmin,r[i].z);
			xmax = MAX(xmax,r[i].x);
			ymax = MAX(ymax,r[i].y);
			zmax = MAX(zmax,r[i].z);
		}
	}
	halo_potent(np,r,penergy);

	for(i=0;i<3;i++){
		ang_mom[i] = 0.L;
	} 
	mrv = 0.L;
	for(i=0;i<np;i++){
		distx = r[i].x - cx; disty = r[i].y - cy; distz = r[i].z - cz;
		distvx = vr[i].vx-cvx; distvy = vr[i].vy-cvy; distvz = vr[i].vz-cvz;
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
	{
		double minE= 1.E20;
		for(i=0;i<np;i++){
			if(penergy[i]+kenergy[i]<minE) {
				minE = penergy[i]+kenergy[i];
				mbp = indx[i];
			}
		}
	}
	tenergy = 0.;
	for(i=0;i<np;i++) tenergy += penergy[i]*0.5 + kenergy[i];
	{
		double a,b;
		a = b = 0.L;
		for(i=0;i<np;i++){
			a += penergy[i];
			b += kenergy[i];
		}
		a = a*0.5L;
		*KW = 2.*b/a+1;
	}
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
		triaxialshape_(&np,x,y,z,&q,&s,&d[0],&v[0][0],&cx,&cy,&cz);
		Free(x);Free(y);Free(z);
	}
	Free(rv);Free(penergy);Free(kenergy);
	return mbp;
}
void getparameter(FILE *wp){
	float bias,rsmooth;
	int npow0;
	float z_init,delta_a;
	float SPHERE_RADIUS,particle_radius,rtheta;
	fscanf(wp,"%f %f\n",&size,&hubble);
	fscanf(wp,"%d %f %f %f %f\n",&npow0,&omep,&omeplam,&bias,&rsmooth);
	fscanf(wp,"%g %d\n",&rng, &nspace);
	fscanf(wp,"%f %f %f\n",&SPHERE_RADIUS,&particle_radius,&rtheta);
	fscanf(wp,"%f %f %f\n",&z_init,&delta_a,&a);
	amax = 1.+z_init;
	redshift = amax/a - 1.;
}


INDXTYPE nhp,maxnhp=30000000;
int sortHaloProp(const void *a, const void *b){
	HaloProp *aa = (HaloProp*)a;
	HaloProp *bb = (HaloProp*)b;
	if(aa->hid < bb->hid ) return -1;
	else if(aa->hid > bb->hid) return 1;
	else return 0;
}
/*
int MAIN_(int argc, char **argv){
*/
int main(int argc, char **argv){
	FILE *fp,*wp,*hid,*mergerwp,*wpbin;
	PtlPos *r;
	PtlVel *vr;
    float peak;
	float *x,*y,*z,*vx,*vy,*vz;
	int tnend,nend,i,j,k;
	long ncount;
	INDXTYPE *indx;
	int id;
    int mpeak,np;
    float xmin,ymin,xmax,ymax,width;
	int haloindex;
	int index,peakindex;
	char infile[100],headerfile[100];
	char wfile[100],mergerfile[100],wbinfile[100];
	int START,inid;
	long nsize;
	int halomergerindex;
	int *nowindex,*beforeindex;
	double kw;

	float tt1,tt2,tt3,tt4;
	int nindex,nmpeak,npeakindex;
	int nhalomergerindex,mtmp;
	int nnowhaloindex,nbeforehaloindex;
	int beforehalomergerindex;
	Halo halo;
	basicparticletype *bp;
	INDXTYPE mx,my,mz,mxmy;
	int finish=-99999;
	int motherrank = 0;
	MPI_Request request;


	HaloProp *hp;


    (void )MPI_Init(&argc,&argv); 
	(void )MPI_Comm_size(MPI_COMM_WORLD,&nid); 
	(void )MPI_Comm_rank(MPI_COMM_WORLD,&myid);


	if(argc != 2) {
		fprintf(stderr,"Error input : FoFhaloQ.exe stepnum \n");
		exit(99);
	}
	(void) Make_Total_Memory(4000L);
	i_potent_spline();
	/*
	den = (float *) Malloc(256*256*256*sizeof(int),PPTR(den));
	*/
	inid = atoi(argv[1]);
	/*
	peak = atof(argv[2]);
	*/
	if(myid==motherrank){
		hp = (HaloProp*)Malloc(sizeof(HaloProp)*maxnhp,PPTR(hp));
		nhp = 0;
		nnowhaloindex = 0;
		sprintf(wfile,"cplFoFHaloQuantities%.5d.dat",inid);
		sprintf(wbinfile,"cplFoFHaloQuantities%.5d.bin",inid);
		wpbin = fopen(wbinfile,"w");
		sprintf(infile,"cplFoF_member_particle.%.5d",inid);
		sprintf(headerfile,"cplFoF_halo_cat.%.5d",inid);
		if((fp = fopen(infile,"r")) == NULL||(hid=fopen(headerfile,"r"))==NULL){
			fprintf(stderr,"error opening file %s && %s\n",infile,headerfile);
			exit(-999);
		}
		fread(&size,sizeof(float),1,hid);
		fread(&hubble,sizeof(float),1,hid);
		fread(&npower,sizeof(float),1,hid);
		fread(&omep,sizeof(float),1,hid);
		fread(&omepb,sizeof(float),1,hid);
		fread(&omeplam,sizeof(float),1,hid);
		fread(&wlam0,sizeof(float),1,hid);
		fread(&wlam1,sizeof(float),1,hid);
		fread(&bias,sizeof(float),1,hid);
		fread(&nx,sizeof(int),1,hid);
		fread(&nspace,sizeof(int),1,hid);
		fread(&amax,sizeof(float),1,hid);
		fread(&astep,sizeof(float),1,hid);
		fread(&anow,sizeof(float),1,hid);
	}
	MPI_Bcast(&size,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&hubble,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&npower,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omep,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omepb,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&omeplam,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&wlam0,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&wlam1,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&bias,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&nx,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&nspace,1,MPI_INT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&amax,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&astep,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	MPI_Bcast(&anow,1,MPI_FLOAT,motherrank,MPI_COMM_WORLD);
	a = anow;


	mx = my = mz =  nx;
	mxmy = mx*my;

	redshift = amax/anow-1;
	rng = nx;
	if(myid==0){
		printf("Now z=%g and rng=%g\n",redshift,rng);
		printf("Now boxsize=%g and hubble=%g\n",size,hubble);
		printf("Now omep=%g omepb=%g omeplam=%g\n",omep,omepb,omeplam);
	}

	nsize = rng*rng*rng/nspace/nspace/nspace;
	INDXTYPE shalonum, rhalonum;



	omepk = 1.-omep-omeplam;

	/*j
   	Hsub = sqrt(omep*pow3(amax/a)+omeplam+omepk*pow2(amax/a));
	*/
	float hofez_(float *, float *, float *, float *, float*, float *);
	Hsub = hofez_(&amax, &anow, &omep, &omeplam, &wlam0, &wlam1);



   	r1kineticfact = size/hubble/rng/amax*a*100.E5*hubble*Hsub;
    r2kineticfact = size/hubble/amax*a*a*100.E5*hubble*Hsub;
	vscale = r2kineticfact/1.E5;
   	sim2real = size/hubble/rng/amax*a*1.E6*pc;
    real2sim = 1./sim2real;
   	pntmass = 3./8./pi/G*sqr(100.E5*hubble)*1.E6*pc;
    pntmass = pntmass*pow3(size/hubble)/pow3(rng/nspace)*omep;
	potentfact = G*pntmass/sim2real;
	index = 0;
	haloindex = 0;
	pscale = anow/amax *size/rng;
	nhp = 0;
	if(myid==motherrank) printf("size=%g rng=%g\n",size,rng);
	MPI_Status mstatus,cstatus;
	int src,dest;
	if(myid==-1){
		int kkk = 1;
		while(kkk){
			kkk = 1;
		}
	}

	if(myid==motherrank){
		INDXTYPE ijstart = 0;
		shalonum = rhalonum = 0;
   		while(fread(&halo,sizeof(Halo),1,hid) == 1){
			double mass;
			double cx,cy,cz;
			if(nhp%10000L ==0) printf("now %ld is being processed\n",nhp);
			{
				i = 0;
				nend = halo.np;
				bp = (basicparticletype *)Malloc(sizeof(basicparticletype)*nend,PPTR(bp));
				fread(bp,sizeof(basicparticletype),nend,fp);

				do {
					MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
					src = mstatus.MPI_SOURCE;
					dest = mstatus.MPI_SOURCE;
					MPI_Recv(&ready,1,MPI_INT,src,READY,MPI_COMM_WORLD,&cstatus);
					if(ready == READY){
						MPI_Send(&nend,sizeof(INDXTYPE),MPI_BYTE,dest,NP_TAG, MPI_COMM_WORLD);
						MPI_Send(bp,nend*sizeof(basicparticletype),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
						MPI_Send(&nhp,sizeof(INDXTYPE),MPI_BYTE,dest,R_TAG,MPI_COMM_WORLD);
						shalonum ++;
					}
					else{
						INDXTYPE inhp;
						MPI_Probe(src,NP_TAG,MPI_COMM_WORLD,&cstatus);
						MPI_Recv(&inhp,sizeof(INDXTYPE),MPI_BYTE,src,NP_TAG,MPI_COMM_WORLD,&cstatus);
#ifdef OLD
						MPI_Recv(hp+inhp,sizeof(HaloProp),MPI_BYTE,src,NP_TAG,MPI_COMM_WORLD,&cstatus);
#else
						MPI_Recv(hp+rhalonum,sizeof(HaloProp),MPI_BYTE,src,NP_TAG,MPI_COMM_WORLD,&cstatus);
#endif
						rhalonum ++;
					}
				} while( ready != READY);
				Free(bp);
				nhp ++;
#ifdef OLD
				if(nhp >= maxnhp){
					maxnhp += 10000000;
					hp = (HaloProp*)Realloc(hp,sizeof(HaloProp)*maxnhp);
				}
#else
				if(rhalonum >= maxnhp){
					maxnhp += 10000000;
					hp = (HaloProp*)Realloc(hp,sizeof(HaloProp)*maxnhp);
				}
#endif
			}
#ifndef OLD
			if(rhalonum%1000000L ==0L) {
				qsort(hp,rhalonum,sizeof(HaloProp),sortHaloProp);
				for(i=0;i<rhalonum;i++){
					if(hp[i].hid != ijstart+i ){
						if(i==0) break;
						ijstart = ijstart+i;
						printf("Writing %d info\n",i);
						fwrite(hp,sizeof(HaloProp),i,wpbin);
						memmove(hp,hp+i,sizeof(HaloProp)*(rhalonum-i));
						rhalonum = rhalonum - i;
						shalonum = shalonum - i;
						break;
					}
				}
			}
#endif
		}
		fclose(hid); fclose(fp);
		printf("complete................... \n");
		fflush(stdout);
		j = 0;
		for(i=rhalonum+1;i<=shalonum;){
			INDXTYPE inhp;
			MPI_Probe(MPI_ANY_SOURCE,READY,MPI_COMM_WORLD,&mstatus);
			MPI_Recv(&ready,1,MPI_INT,mstatus.MPI_SOURCE,READY,
					MPI_COMM_WORLD,&cstatus);
			if(ready == WRITING){
				i++;
				MPI_Probe(mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD,&cstatus);
				MPI_Recv(&inhp,sizeof(INDXTYPE),MPI_BYTE,mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD, &cstatus);
#ifdef OLD
				MPI_Recv(hp+inhp,sizeof(HaloProp),MPI_BYTE,mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD, &cstatus);
#else
				MPI_Recv(hp+rhalonum,sizeof(HaloProp),MPI_BYTE,mstatus.MPI_SOURCE,NP_TAG,MPI_COMM_WORLD, &cstatus);
				rhalonum ++;
#endif
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
		printf("finished gathering data \n");fflush(stdout);
	}
	else {
		ready = READY;
		MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
		MPI_Recv(&nend,sizeof(INDXTYPE),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
		while(nend != finish){
			bp = (basicparticletype*)Malloc(sizeof(basicparticletype)*nend,PPTR(bp));
			MPI_Recv(bp,sizeof(basicparticletype)*nend,MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			MPI_Recv(&nhp,sizeof(INDXTYPE),MPI_BYTE,motherrank,R_TAG,MPI_COMM_WORLD,&cstatus);
			r = (PtlPos *) Malloc(nend*sizeof(PtlPos),PPTR(r));
			vr = (PtlVel *) Malloc(nend*sizeof(PtlVel),PPTR(vr));
			indx = (INDXTYPE *) Malloc(nend*sizeof(INDXTYPE),PPTR(indx));
			for(j=0;j<nend;j++){
				r[j].x = bp[j].x;
				r[j].y = bp[j].y;
				r[j].z = bp[j].z;
				vr[j].vx = bp[j].vx;
				vr[j].vy = bp[j].vy;
				vr[j].vz = bp[j].vz;
				indx[j] = bp[j].indx;
			}
			if(nend == 0) {
				continue;
			}
			cx = cy = cz = 0;
			for(j=0;j<nend;j++){
				/* WARNING */
				cx += r[j].x;
				cy += r[j].y;
				cz += r[j].z;
			}
			cx = cx / (double) nend; 
			cy = cy / (double) nend; 
			cz = cz / (double) nend;
			INDXTYPE mbp = halo_properties_kw(nend,r,vr,indx,&kw);
			{
				float *px,*py,*pz;
				px = (float *) Malloc(sizeof(float)*10000,PPTR(px));
				py = (float *) Malloc(sizeof(float)*10000,PPTR(py));
				pz = (float *) Malloc(sizeof(float)*10000,PPTR(pz));
				r1[0] = v[0][0]; r1[1] = v[0][1]; r1[2] = v[0][2];
				r2[0] = v[1][0]; r2[1] = v[1][1]; r2[2] = v[1][2];
				r3[0] = v[2][0]; r3[1] = v[2][1]; r3[2] = v[2][2];
				radius = find_r(r,nend,q,s,r1,r2,r3)*size/rng;
				Free(pz);
				Free(py);
				Free(px);
			}
			float mass = nend*pntmass/onesolarmass;
			{
				double vvx,vvy,vvz;
				cvx = cvy = cvz = 0.;
				for(j=0;j<nend;j++){
					cvx += vr[j].vx;
					cvy += vr[j].vy;
					cvz += vr[j].vz;
				}
				cvx = cvx/(double)nend;
				cvy = cvy/(double)nend;
				cvz = cvz/(double)nend;
				sigv = 0;
				for(j=0;j<nend;j++){
					vvx = (vr[j].vx-cvx)*r2kineticfact/1.E5 + (r[j].x-cx)*r1kineticfact/1.E5;
					vvy = (vr[j].vy-cvy)*r2kineticfact/1.E5 + (r[j].y-cy)*r1kineticfact/1.E5;
					vvz = (vr[j].vz-cvz)*r2kineticfact/1.E5 + (r[j].z-cz)*r1kineticfact/1.E5;
					sigv += vvx*vvx+vvy*vvy+vvz*vvz;
				}
				sigv = sigv/(double) nend;
				sigv = sqrt(sigv);
			}
			/* pvr : peculiar velocity in km/sec */
			pvr = sqrt(cvx*cvx+cvy*cvy+cvz*cvz)*r2kineticfact/1.E5;
			if(isnan(para_spin) != 0) {
				para_spin = -9999.;
			}
			cvx = cvx*vscale;
			cvy = cvy*vscale;
			cvz = cvz*vscale;
			{
				double ampang;
				ampang = pow2(ang_mom[0])+pow2(ang_mom[1])+pow2(ang_mom[2]);
				ampang = sqrt(ampang);
				ang_mom[0] = ang_mom[0]/ampang;
				ang_mom[1] = ang_mom[1]/ampang;
				ang_mom[2] = ang_mom[2]/ampang;
			}
			Free(indx);Free(vr);Free(r);Free(bp);
			HaloProp ahp;
			ahp.x = cx*size/rng;
			ahp.y = cy*size/rng;
			ahp.z = cz*size/rng;
			ahp.mass = mass;
			ahp.hid = nhp;
			ahp.mbp = mbp;
			ahp.paraspin = para_spin;
			ahp.sigv = sigv;
			ahp.rotation = rotation;
			ahp.radius = radius;
			ahp.q = q;
			ahp.s = s;
			ahp.vx = cvx;
			ahp.vy = cvy;
			ahp.vz = cvz;
			for(j=0;j<3;j++) {
				for(k=0;k<3;k++){
					ahp.v[j][k] = v[j][k];
				}
				ahp.ang_mom[j] = ang_mom[j];
			}

			ready = WRITING;
			MPI_Issend(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD, &request);
			MPI_Wait(&request,&mstatus);
			MPI_Send(&nhp,sizeof(INDXTYPE),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD);
			MPI_Send(&ahp,sizeof(HaloProp),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD);

			ready = READY;
			MPI_Send(&ready,1,MPI_INT,motherrank,READY,MPI_COMM_WORLD);
			MPI_Recv(&nend,sizeof(INDXTYPE),MPI_BYTE,motherrank,NP_TAG,MPI_COMM_WORLD,&cstatus);
		}
	}
	if(myid==motherrank){
		/*
		INDXTYPE ii;
		wp = fopen(wfile,"w");
		HaloProp *pt = hp;
		for(ii=0;ii<nhp;ii++){
			fprintf(wp,"%g %f %f %f %ld %ld %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g  %g %g %g %g %g %g %g\n",
					pt->mass,pt->x,pt->y,pt->y,pt->z,ii,pt->mbp,pt->paraspin,pt->sigv,pt->rotation,
					pt->radius,pt->q,pt->s,pt->v[0][0],pt->v[0][1],pt->v[0][2],
					pt->v[1][0],pt->v[1][1],pt->v[1][2],
					pt->v[2][0],pt->v[2][1],pt->v[2][2],
					pt->ang_mom[0],pt->ang_mom[1],pt->ang_mom[2],pt->vx,pt->vy,pt->vz);
			pt++;
		}
		fflush(wp);
		fclose(wp);
		*/
#ifdef OLD
		qsort(hp,nhp,sizeof(HaloProp),sortHaloProp);
		/*
		fwrite(&nhp,sizeof(INDXTYPE),1,wpbin);
		*/
		fwrite(hp,sizeof(HaloProp),nhp,wpbin);
#else
		qsort(hp,rhalonum,sizeof(HaloProp),sortHaloProp);
		fwrite(hp,sizeof(HaloProp),rhalonum,wpbin);
#endif
		fclose(wpbin);
	}
	Free(den);
	MPI_Finalize();
	return 0;
}
float piinhin3;
float h;
void findden(int np, PtlPos *r, int nxwidth,int nywidth,int nzwidth){
	int i,j,k;
	int nbuff;
	float tmp2h;
	
	piinhin3 = 1./3.141592653;
	h = 1.;
	nbuff = 2*h+1;
	tmp2h = 2.*h;
        for(i=0;i<np;i++){
                float xi,yi,zi,xf,yf,zf;
                int mxi,myi,mzi;
                int mxf,myf,mzf;
                float wr,wresult,gx,gy,gz;
                int ii,ij,ik;
                int iii,iij,iik;
                float xp,yp,zp;
                xp = r[i].x;
                yp = r[i].y;
                zp = r[i].z;
                xi = xp - tmp2h; yi = yp - tmp2h; zi = zp - tmp2h;
                xf = xp + tmp2h; yf = yp + tmp2h; zf = zp + tmp2h;
                mxi = rint(xi) + 1; myi = rint(yi) + 1; mzi = rint(zi) + 1;
                mxf = rint(xf); myf = rint(yf); mzf = rint(zf);
                for(ik = mzi;ik<=mzf;ik++)
                for(ij = myi;ij<=myf;ij++)
                for(ii = mxi;ii<=mxf;ii++) {
                        float tmpx,tmpy,tmpz;
                        gx = ii - 0.5; gy = ij - 0.5; gz = ik - 0.5;
                        tmpx = gx-xp;tmpx = tmpx*tmpx;
                        tmpy = gy-yp;tmpy = tmpy*tmpy;
                        tmpz = gz-zp;tmpz = tmpz*tmpz;
                        wr = sqrt(tmpx + tmpy + tmpz);
                        /*
                        wr = w_(&wr);
                        */
                        wresult = wc(wr);
                        den[ii+ij*nxwidth+ik*nxwidth*nywidth] += wresult;
                }
        }
}
float wc(float r){
        float w,x;
        x = r/h;
        if(x < 1.)
                w = piinhin3*(1.-1.5*x*x+0.75*x*x*x);
        else if(x < 2.)
                w = piinhin3*0.25*(2.-x)*(2.-x)*(2.-x);
        else w = 0;
        return w;
}
