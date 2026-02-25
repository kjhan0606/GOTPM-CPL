/* My tree correction algorithm */
/* Modified from size_t to long integer due to the problem
 * occuring when comparison is made using signed and unsigned integers
 * 10/12/2005
 *
 * OpenMP directives are inserted to allow the parallel computing on SMP.
 * Due to the canonical form of for loop parallelization, I have to
 * use integer index rather than long index.
 * But long index will be promising for the future implementation
 * to allow very large simulation particles 
 * 09/03/2006
 *
 * Arbitrary zwith of local domain is now allowed 13/03/2006: being tested 16/03/2006
 * */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<omp.h>
#include "mpi.h"
#include "Memory.h"
#include "pmheader.h"
#include "kjhtree.h"
#include "force_spline.h"

#ifdef USE_GPU
#include "nvgpu.h"
void gputreeforce(particle *,int, TStruct *,int , TPtlStruct *,int, float, int, int,int);
void Direct_Nbody(particle *,int, TPtlStruct *,int , int, int,int);
#endif

/*
float cputime0[10],cputime1[10];
float gettime();
*/
#define min(a,b) ((a)<(b) ? (a):(b))
#define max(a,b) ((a)>(b) ? (a):(b))
#define CellWidth 4.L /* The cellwidth is preferred to be same as rsphere */
#define invCellWidth 0.25L /* This is 1/CellWidth . */
#define rsphere 4.
#define nstride 1
float rspheresq;
float maintreetime;
typedef struct ComGhost{
#ifdef MULTIMASS
	float mass;
#endif
#ifdef XYZDBL
	double x,y,z;
#else
	float x,y,z;
#endif
	struct ComGhost *next;
} ComGhost;

static int initflag=0;
long mx,my,mz;
treeparticletype **BasicCell;
long *NBasicCell;
int nullfct0() {
return 0;
}
int nullfct1() {
return 1;
}
#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif

/*
long npos,ngrv;
particle *pos,*tpos;
TPtlStruct *grv,*tgrv;
*/

void Ghost(long,float,float,float,int,int,int,int,float,float);
void i_force_spline(long,float);
int dest,src;
int idest,isrc;
int myid,nid;
float theta2;
void ParticleCellLinkedList(treeparticletype *particles, long np, float zstart,
		int nx,int ny,int nz,int nspace){
	long i,ix,iy,iz,mpixel;
	long mxm1,mym1,mzm1;
	treeparticletype *p,*tp;
	for(i=0;i<mx*my*mz;i++){
		BasicCell[i] = NULL;
		NBasicCell[i] = 0;
	}

	mxm1 = mx - 1;
	mym1 = my - 1;
	mzm1 = mz - 1;
	tp = particles;
	for(i=0;i<np;i++){
		ix = (long)(XofP(tp)*invCellWidth);
		iy = (long)(YofP(tp)*invCellWidth);
		iz = (long)((ZofP(tp)-zstart)*invCellWidth);
		/* This is for the case of rounding-off issue */
		ix = min(ix,mxm1); iy = min(iy,mym1); iz = min(iz,mzm1);
		mpixel = ix+mx*(iy+my*iz);
		NBasicCell[mpixel] ++;
		p = BasicCell[mpixel];
		BasicCell[mpixel] = tp;
		tp->next = p;
		tp++;
	}
}
int mbx;
void BoundaryCellLinkedList(ComGhost *bparticles, long nbp, 
		ComGhost **BoundaryCell, long *NBoundaryCell){
	long i,ix,iy,iz,mpixel;
	ComGhost *p,*tp;
	for(i=0;i<mbx;i++){
		BoundaryCell[i] = NULL;
		NBoundaryCell[i] = 0;
	}
	tp = bparticles;
	for(i=0;i<nbp;i++){
		ix = (long)(tp->x/CellWidth);
		p = BoundaryCell[ix];
		BoundaryCell[ix] = tp;
		tp->next = p;
		NBoundaryCell[ix] ++;
		tp++;
	}
}
Box CheckBoundingBox(long ngrv,TPtlStruct *grv){
	float xmin,xmax,ymin,ymax,zmin,zmax;
	float width;
	long i;
	int si;
	Box box;
	xmin = ymin = zmin = 1.e25;
	xmax = ymax = zmax = -1.e25;
	for(i=0;i<ngrv;i++){
		xmin = min(xmin,grv[i].r[0]);
		ymin = min(ymin,grv[i].r[1]);
		zmin = min(zmin,grv[i].r[2]);
		xmax = max(xmax,grv[i].r[0]);
		ymax = max(ymax,grv[i].r[1]);
		zmax = max(zmax,grv[i].r[2]);
	}
	width = xmax-xmin;
	width = max(width,ymax-ymin);
	width = max(width,zmax-zmin);
	if(width !=0){
		box.x = xmin-width*1.e-6;
		box.y = ymin-width*1.e-6;
		box.z = zmin-width*1.e-6;
		box.width = width*(1.+2.e-6);
	}
	else {
		box.x = xmin-1.e-5;
		box.y = ymin-1.e-5;
		box.z = zmin-1.e-5;
		box.width = 2.e-5;
	}
	return box;
}
/*
void WringGrvInfo(long ncell){
	treeparticletype *p;
	long ntmp;
	p = BasicCell[ncell];
	ntmp = 0;
	while(p){
		grv[ntmp].r[0] = p->x;
		grv[ntmp].r[1] = p->y;
		grv[ntmp].r[2] = p->z;
		p = p->next;
		ntmp ++;
	}
}
*/
/* Dump2Grv is designed to copy particle data in grv format considering isolated boundaries */
long GrvCountParticle(long ix0,long iy0,long iz0,long ixm,
		long ixp,long iym,long iyp,long izm,long izp) {
	long ix,iy,iz;
	long ixi,iyi,izi,ixf,iyf,izf;
	long mcell,ncount;
	long zshift,yzshift,mxy;
	/* for isolated boundary conditions */
	ixi = max(ix0+ixm,0);
	iyi = max(iy0+iym,0);
	izi = max(iz0+izm,0);
	ixf = min(ix0+ixp,mx);
	iyf = min(iy0+iyp,my);
	izf = min(iz0+izp,mz);

	ncount = 0;
	mxy = mx*my;
	for(iz=izi;iz<izf;iz++){
		zshift = mxy*iz;
		for(iy=iyi;iy<iyf;iy++){
			yzshift = mx*iy+zshift;
			for(ix=ixi;ix<ixf;ix++){
				mcell = ix+yzshift;
				ncount += NBasicCell[mcell];
			}
		}
	}
	return ncount;

}
long Dump2Grv(long ix0,long iy0,long iz0,long ixm,long ixp,
		long iym,long iyp,long izm,long izp,TPtlStruct *grv, 
		int nx,int ny,int nz,int nspace, double xtran, 
		double ytran, double ztran) {
	long ix,iy,iz,mcell;
	long ixi,iyi,izi,ixf,iyf,izf;
	long ntmp,ngrv;
	long zshift,yzshift,mxy;
	TPtlStruct *tgrv;
	treeparticletype *p;
	tgrv = grv;
	ngrv = 0;
	/* for isolated boundary conditions */
	ixi = max(ix0+ixm,0);
	iyi = max(iy0+iym,0);
	izi = max(iz0+izm,0);
	ixf = min(ix0+ixp,mx);
	iyf = min(iy0+iyp,my);
	izf = min(iz0+izp,mz);

	mxy = mx*my;
	for(iz=izi;iz<izf;iz++){
		zshift = mxy*iz;
		for(iy=iyi;iy<iyf;iy++){
			yzshift = mx*iy+zshift;
			for(ix=ixi;ix<ixf;ix++){
				mcell = ix+yzshift;
				p=BasicCell[mcell];
				while(p){
					grv[ngrv].type = TYPE_PTL;
					grv[ngrv].r[0] = XofP(p)+xtran;
					grv[ngrv].r[1] = YofP(p)+ytran;
					grv[ngrv].r[2] = ZofP(p)+ztran;
#ifdef MULTIMASS
					grv[ngrv].mass = p->mass;
#else
					grv[ngrv].mass = 1.;
#endif
					p = p->next;
					ngrv++;
				}
			}
		}
	}
	return ngrv;
}
long PosCountParticle(long ix0,long iy0,long iz0,long imx,
		long ipx,long imy,long ipy,long imz,long ipz) {
	long ncount,mcell;
	long ix,iy,iz;
	long ix2,iy2,iz2;
	long zshift,yzshift,mxy;
	ncount = 0;
	mxy = mx*my;
	for(iz=iz0+imz;iz<iz0+ipz;iz++){
		if(iz < 0 || iz >= mz) continue;
		else iz2 = iz;
		zshift = mxy*iz2;
		for(iy=iy0+imy;iy<iy0+ipy;iy++){
			iy2 = (iy+my)%my;
			yzshift = mx*iy2+zshift;
			for(ix=ix0+imx;ix<ix0+ipx;ix++){
				ix2 = (ix+mx)%mx;
				mcell = ix2+yzshift;
				ncount += NBasicCell[mcell];
			}
		}
	}
	return ncount;
}

long Dump2PeriodicPosExceptZ(long ix0,long iy0,long iz0,long ixm,
		long ixp, long iym, long iyp, long izm, long izp,
		int nx,int ny,int nz,int nspace, particle *pos,
		double xtran, double ytran, double ztran){
	long ix,iy,iz;
	long ix2,iy2,iz2;
	long mcell;
	long zshift,yzshift,mxy;
	double lpx,lpy;
	long npos;
	treeparticletype *p;
	npos = 0;
	mxy = mx*my;
	for(iz=iz0+izm;iz<iz0+izp;iz++){
		if(iz < 0 || iz >= mz) continue;
		else iz2 = iz;
		zshift = mxy*iz2;
		for(iy=iy0+iym;iy<iy0+iyp;iy++){
			iy2 = (iy+my)%my;
			if(iy > iy2) lpy = ny;
			else if(iy < iy2) lpy = -ny;
			else lpy = 0.;
			yzshift = mx*iy2+zshift;
			for(ix=ix0+ixm;ix<ix0+ixp;ix++){
				ix2 = (ix+mx)%mx;
				mcell = ix2+yzshift;
				p=BasicCell[mcell];
				if(ix > ix2) lpx = nx;
				else if(ix < ix2) lpx = -nx;
				else lpx = 0.;
				while(p){
#ifdef INDTIME
					if(isnowstep(p->tflag,maxsubT))
#endif
					{
						pos[npos].x = XofP(p) + lpx + xtran;
						pos[npos].y = YofP(p) + lpy + ytran;
						pos[npos].z = ZofP(p) + ztran;
						pos[npos].bp = p;
					}
					npos ++;
					p = p->next;
				}
			}
		}
	}
	return npos;
}
/*
void determine_mpi_long(){
	if(sizeof(long)==sizeof(unsigned int)) {
		MPI_long = MPI_UNSIGNED;
	}
	else if(sizeof(long) == sizeof(unsigned long)){
		MPI_long = MPI_UNSIGNED_LONG;
	}
}
long MaxAllowableCell(){
	long nfree,ncell,needed,additional,nscale;

	nfree = freespace();
	needed = (sizeof(TStruct)+sizeof(TPtlStruct)+sizeof(particle)*27)*npmaxcell;

	return ncell;
}
*/
/* v(t+dt/2) = fact1 * v(t-dt/2) + fact2 + a(t-dt/2) 
 * r(t+dt) = r(t-dt) + v(t+dt/2)*pfact
 */
void UpdateVelocityUsingForceByParticles(long np, float vfact2,particle *pos){
	long i;
#ifdef _OPENMP
	int si;
#pragma omp parallel for 
	for(si=0;si<(int)np;si++){
		(pos[si].bp)->vx += vfact2*pos[si].ax*simpar.pcorr;
		(pos[si].bp)->vy += vfact2*pos[si].ay*simpar.pcorr;
		(pos[si].bp)->vz += vfact2*pos[si].az*simpar.pcorr;
	}
#else
	for(i=0;i<np;i++){
		(pos[i].bp)->vx += vfact2*pos[i].ax*simpar.pcorr;
		(pos[i].bp)->vy += vfact2*pos[i].ay*simpar.pcorr;
		(pos[i].bp)->vz += vfact2*pos[i].az*simpar.pcorr;
	}
#endif
}
void DirectSummation(long ngrv,TPtlStruct *grv,long npos,particle *pos){
	float tmpx,tmpy,tmpz,fplmf;
	float px,py,pz,dist2,ptlmass,dist;
	long i,j,ntmp;
	for(i=0;i<npos;i++)
		pos[i].ax = pos[i].ay = pos[i].az = 0.;
	
	for(i=0;i<npos;i++){
		px = pos[i].x; py = pos[i].y; pz = pos[i].z;
		for(j=0;j<ngrv;j++){
			tmpx = px - grv[j].r[0];
			tmpy = py - grv[j].r[1];
			tmpz = pz - grv[j].r[2];
			dist2 = tmpx*tmpx + tmpy*tmpy + tmpz*tmpz;
			if(dist2 <= rspheresq){
				dist = SQRT(dist2);
				ntmp = dist*invran2nran;
				fplmf = forcecorrectdiff(ntmp,0) + forcecorrectslope(ntmp,0) * (dist-ntmp*ran2nran);  
				pos[i].ax += tmpx* fplmf;
				pos[i].ay += tmpy* fplmf;
				pos[i].az += tmpz* fplmf;
			}
		}
	}
}
void PartialUpdateAccelUsingForceByParticles(float vfact2,long npos,particle *pos,float mass){
	long i;
	treeparticletype *bp;
	float xmin,ymin,zmin;
	float xmax,ymax,zmax;
	for(i=0;i<npos;i++){
		bp = pos[i].bp;
		bp->vx += vfact2*pos[i].ax*mass*simpar.pcorr;
		bp->vy += vfact2*pos[i].ay*mass*simpar.pcorr;
		bp->vz += vfact2*pos[i].az*mass*simpar.pcorr;
	}
}
/*
void AllUpdateVelocityAndPosition(treeparticletype *particles,long np,
		int nx, int ny, int nz,
		float pfact, float vfact1, float vfact2){
	long i;
	treeparticletype *tp;
	tp = particles;
	for(i=0;i<np;i++){
		tp->vx = vfact1* tp->vx + vfact2 * tp->ax;
		tp->vy = vfact1* tp->vy + vfact2 * tp->ay;
		tp->vz = vfact1* tp->vz + vfact2 * tp->az;
		tp->x += pfact * tp->vx;
		tp->y += pfact * tp->vy;
		tp->z += pfact * tp->vz;
		tp++;
	}
	tp = particles;
	for(i=0;i<np;i++){
		tp->x = fmod(tp->x+nx,nx);
		tp->y = fmod(tp->y+ny,ny);
		tp->z = fmod(tp->z+nz,nz);
	}
}
*/

long npmaxcell;
float treecorrection(treeparticletype *particles,long np){
	float mass,zstart,zheight;
	int nx,ny,nz,nspace;
	float vfact1,vfact2,pfact,theta;
	long i,j,k;
	long ix0,iy0,iz0;
	long ix2,iy2,iz2;
	long jxp,jyp,jzp;
	long ncell,mcell,mp;
	long npos,ngrv;
	long mpos,mgrv;
	double xtran,ytran,ztran;
	particle *pos;
	TPtlStruct *grv;
	TStruct *TREECELL;
	treeparticletype *p;
	Box box;
#ifdef _OPENMP
	int si; /* for integer*4 value of loop index in OpenMP*/
#endif
	mass = simpar.nspace*simpar.nspace*simpar.nspace;
	zstart = simpar.zmin;
	zheight = (simpar.zmax-simpar.zmin);
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	vfact1 = simpar.fact1;
	vfact2 = simpar.fact2;
	pfact = simpar.pfact;
	theta = simpar.theta;

	rspheresq = rsphere*rsphere;
	/*
	determine_mpi_long();
	*/
	MPI_Comm_size(MPI_COMM_WORLD,&nid);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);


	theta2 = theta*theta;



	mx = (long)ceil((double)nx/(double)CellWidth);
	my = (long)ceil((double)ny/(double)CellWidth);
	mz = (long)ceil((double)zheight/(double)CellWidth);

	if(initflag==0) {
		if(nx != ny && nx!=nz) {
			if(myid==0) {
				fprintf(stderr,"Error in TreeCorrection\n");
				fprintf(stderr,"since nx != ny, we don't know how to set force scaling \n");
			}
			exit(99);
		}
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid==0) printf("#######################################\n");
		if(myid==0) printf("P%d Now initializing force array\n",myid);
		if(myid==0) printf("#######################################\n\n");
		MPI_Barrier(MPI_COMM_WORLD);
		i_force_spline(nx,rsphere);
		initflag = 1;
	}

	p = particles;
	BasicCell = (treeparticletype **) 
		Malloc(sizeof(treeparticletype *)*mx*my*mz,PPTR(BasicCell));
	NBasicCell = (long *) Malloc(sizeof(long)*mx*my*mz,PPTR(NBasicCell));
	ParticleCellLinkedList(particles,np,zstart,nx,ny,nz,nspace);
#ifdef DEBUG 
	printf("P%d has np=%ld mx=%ld my=%ld mz=%ld memory for basicell = %ld zstart=%g vfact2=%g\n",
			myid,np, mx,my,mz, sizeof(treeparticletype *)*mx*my*mz,zstart,vfact2);
#endif
	{
		long nposmax,ngrvmax;
		nposmax = 0;
		for(i=0;i<mx*my*mz;i++) nposmax = max(nposmax,NBasicCell[i]);
		npmaxcell = nposmax;
		mpos = nposmax * (nstride+2)*(nstride+2)*(nstride+2);
		mgrv = nposmax * nstride*nstride*nstride;

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			if(omp_get_thread_num()==0) {
				pos = (particle *)Malloc(sizeof(particle)*mpos*omp_get_num_threads(),PPTR(pos));
				grv = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*mgrv*omp_get_num_threads(),PPTR(grv));
				TREECELL = (TStruct *)Malloc(sizeof(TStruct)*mgrv*omp_get_num_threads(),PPTR(TREECELL));
			}
		}
	}
#ifdef USE_GPU
	GPUDIRECTSUM = 2000 - 5000*(simpar.theta-0.2);
	GPUDIRECTSUM = min(GPUDIRECTSUM, 3000);
	GPUDIRECTSUM = max(GPUDIRECTSUM,  800);
#endif



	for(iz0=0;iz0<mz;iz0+=nstride){
		jzp = min(iz0+nstride,mz)-iz0+1;
		/* This OpenMP struct could have the data racing conditions 
		 * if "chunk" = 1 and 2 or possibly 3 (?).
		 * And "my" should be divisible by "nthreads" and */
#ifdef _OPENMP
#pragma omp parallel private(iy0,jyp,jxp,ix0,ngrv,npos,i,box,xtran,ytran,ztran)
#endif
		{
			int siy0;
			int idthread,nthreads;
			long myiy0start,myiy0end,chunk;
			int lcount=0;
			particle *threadpos;
			TPtlStruct *threadgrv;
			TStruct *threadTREECELL;
			idthread = omp_get_thread_num();
			nthreads = omp_get_num_threads();
			threadpos = pos + mpos*idthread;
			threadgrv = grv + mgrv*idthread;
			threadTREECELL = TREECELL + mgrv*idthread;
			chunk = (my+nthreads-1)/nthreads;
			myiy0start = chunk*idthread;
			myiy0end = min(myiy0start+chunk,my);
			ztran = -(iz0*CellWidth+zstart);
			/*
			if(myid==0) printf("P%d has %ld %ld for %d\n",myid,myiy0start,myiy0end,idthread);
			*/
#ifdef _OPENMP
			for(iy0=myiy0start;iy0<myiy0end;iy0+=nstride){
#else
			for(iy0=0;iy0<my;iy0+=nstride){
#endif
				ytran = -iy0*CellWidth;
				jyp = min(iy0+nstride,my)-iy0+1;
				for(ix0=0;ix0<mx;ix0+=nstride){
					jxp = min(ix0+nstride,mx)-ix0+1;
					xtran = -ix0*CellWidth;
					/* Making tree buffer */
					ngrv = Dump2Grv(ix0,iy0,iz0,0,nstride,0,nstride,0,nstride,threadgrv,
							nx,ny,nz,nspace, xtran,ytran,ztran);
					if(ngrv > 0){
						npos = Dump2PeriodicPosExceptZ(ix0,iy0,iz0,-1,jxp,-1,jyp,-1,jzp,
								nx,ny,nz,nspace,threadpos,
								xtran,ytran,ztran);
#ifndef USE_GPU
						if(ngrv >= DIRECTSUM){
							box.x = box.y = box.z = 0;
							box.width = nstride*CellWidth;
							Make_Tree(threadTREECELL,threadgrv,ngrv,box,theta);
							for(i=0;i<npos;i++)
								treeforce(threadpos+i,theta2,threadTREECELL,threadgrv,rspheresq);
						}
						else DirectSummation(ngrv,threadgrv,npos,threadpos);
#else
						/*
						TIMER_START(45);
						*/
						if(ngrv <= GPUCPUDIRECTSUM){
							DirectSummation(ngrv,threadgrv,npos,threadpos);
						}
						else if(ngrv <= GPUDIRECTSUM){
							Direct_Nbody(threadpos,npos, threadgrv,ngrv , GPUSPERNODE, myid,nx);
						}
						else { /* This is GPU version of tree sum */
							int ncell;
							box.x = box.y = box.z = 0;
							box.width = nstride*CellWidth;
							ncell = Make_Tree(threadTREECELL,threadgrv,ngrv,box,theta);
							gputreeforce(threadpos,npos,threadTREECELL,ncell,threadgrv,ngrv,
								rspheresq,GPUSPERNODE,myid,nx);
						}
						/*
						TIMER_STOP(45);
						if(lcount ==1) {
							printf("P%d Total GPU Time = %g for ngrv= %d\n",myid,ELAPSED_TIME(45),ngrv);
							MPI_Finalize();
							exit(99);
						}
						lcount ++;
						*/
#endif
						PartialUpdateAccelUsingForceByParticles(vfact2,npos,threadpos,mass);
					}
				}
			}
		}
	}
#ifdef _OPENMP
#pragma omp parallel
#endif
	{
			if(omp_get_thread_num()==0) {
					Free(TREECELL);
					Free(grv);
					Free(pos);
			}
	}


	TIMER_START(37);
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid==0) printf("P%d now entering ghost routine 1\n",myid);
		Ghost(1,zstart,zheight,vfact2,nx,ny,nz,nspace,mass,theta);
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid==0) printf("P%d now entering ghost routine -1\n",myid);
		Ghost(-1,zstart,zheight,vfact2,nx,ny,nz,nspace,mass,theta);
		MPI_Barrier(MPI_COMM_WORLD);
	}
	TIMER_STOP(37);
    if(myid==0) fprintf(stdout,"Ghost Tree. CPU= %f \n", ELAPSED_TIME(37));

	Free(NBasicCell);
	Free(BasicCell);
	return maintreetime;
}

/* if iflag == 1 --> upper sending 
 * else if iflag == -1 --> down sending
 * else --> Error
 * */
void Ghost(long iflag,float zstart,float zheight,float vfact2,
		int nx, int ny, int nz,int nspace,float mass, float theta){
	long i,j,k;
	long nrghost,nsghost,mcell;
	long ix,iy,iz;
	long ix0,iy0,izg0;
	long jxp,jyp,jzm,jzp;
	long ix1,iy1,izg1;
	long iz0pos;
	long ncount,ngrv,npos;
	particle *pos;
	TPtlStruct *grv;
	int sndtag,recvtag;
	MPI_Status status;
	float z_start,z_stop;
	float lpz;
	treeparticletype *p;
	TStruct *TREECELL;
	Box box;
	ComGhost *sghost,*rghost;
	ComGhost **BoundaryCell;
	long *NBoundaryCell;
	long mpos,mgrv;
	double xtran,ytran,ztran;
#ifdef _OPENMP
	int si;
#endif
	float *zw;
	float twidth;
	long nneighbor,tnneighbor; /* starts from 1 */
	long nget, nput;
	float zminput,zmaxput,zwsum;
	long zshift,yzshift,mxy;
	int *displs;


	sndtag = recvtag = 1;
	if(iflag == 1) {
		dest = (myid+1+nid)%nid;
		src = (myid-1+nid)%nid;
		izg1 = mz-1;/* plus slice (being included) for gravitation data */
		if(fmod(zheight,rsphere)==0.L) izg0 = izg1;
		else izg0 = izg1-1;
		if(mz==1) izg0 = 0;
		z_start = zstart+zheight-rsphere;
		z_stop = zstart+zheight;
		iz0pos = 0;
		jzm = 0; /* minus slice for position data */
		jzp = 1; /* plus slice for position data */
	}
	else {
		dest = (myid-1+nid)%nid;
		src = (myid+1+nid)%nid;
		izg1 = izg0 = 0;
		z_start = zstart;
		z_stop = zstart+rsphere;
		iz0pos = mz-1;
		if(fmod(zheight,rsphere)==0.L){
			jzm = 0;
			jzp = 1;
		}
		else{
			jzm = -1;
			jzp = 1;
		}
	}

	if(iflag ==1 && myid==nid-1) lpz = -nz;
	else if(iflag == -1 && myid== 0) lpz = nz;
	else lpz = 0;
	{ /* Check whether communication beyond contact neighbors is needed.
		 It outputs the number of communications (nneighbor,
		 and its globale maxima tnneighbor) beyond boundary */
		zw=(float*)Malloc(sizeof(float)*nid,PPTR(zw));

		if(myid==0) {
			zw[0] = zheight;
			for(i=1;i<nid;i++)
				MPI_Recv(zw+i,1,MPI_FLOAT,i,i,MPI_COMM_WORLD,&status);
		}
		else 
			MPI_Send(&zheight,1,MPI_FLOAT,0,myid,MPI_COMM_WORLD);

		MPI_Bcast(zw,nid,MPI_FLOAT,0,MPI_COMM_WORLD);
		nneighbor = 1; /* starting from 1 */
		twidth = 0;
		for(i=1;i<nid;i++){
			twidth += zw[(myid+i*iflag+nid)%nid];
			if(twidth < rsphere) nneighbor ++;
			else break;
		}
		MPI_Reduce(&nneighbor,&tnneighbor,1,MPI_LONG,MPI_MAX,0,MPI_COMM_WORLD);
		MPI_Bcast(&tnneighbor,1,MPI_LONG,0,MPI_COMM_WORLD);
	}
	mxy = mx*my;
	ztran = -(iz0pos*CellWidth+zstart);
	for(iy0=0;iy0<my;iy0+=nstride){
		iy1 = min(my,iy0+nstride);
		jyp = iy1-iy0+1;
		sghost = (ComGhost*)Malloc(CheckAvailableMemory()-100,PPTR(sghost));
		ncount = 0;
		ytran = -iy0*CellWidth;
		for(iz=izg0;iz<=izg1;iz++){
			zshift = mx*my*iz;
			for(iy=iy0;iy<iy1;iy++){
				yzshift = mx*iy+zshift;
				for(ix=0;ix<mx;ix++){
					mcell = ix + yzshift;
					p = BasicCell[mcell];
					while(p){
#ifdef XYZDBL
						double zzz;
						zzz = ZofP(p);
						if(zzz >=z_start && zzz < z_stop){
#ifdef MULTIMASS
							sghost[ncount].mass = p->mass;
#endif
							sghost[ncount].x = XofP(p);
							sghost[ncount].y = YofP(p);
							sghost[ncount].z = zzz;
							ncount++;
						}
#else
						if(p->z >=z_start && p->z < z_stop){
#ifdef MULTIMASS
							sghost[ncount].mass = p->mass;
#endif
							sghost[ncount].x = p->x;
							sghost[ncount].y = p->y;
							sghost[ncount].z = p->z;
							ncount++;
						}
#endif
						p= p->next;
					}
				}
			}
		}
		nsghost = ncount;
		sghost = (ComGhost*)Realloc(sghost,sizeof(ComGhost)*nsghost);
		MPI_Sendrecv(&nsghost,1,MPI_LONG,dest,sndtag,
				&nrghost,1,MPI_LONG,src,recvtag,MPI_COMM_WORLD,
				&status);
		rghost = (ComGhost*)Malloc(sizeof(ComGhost)*nrghost,PPTR(rghost));
		MPI_Sendrecv(sghost,sizeof(ComGhost)*nsghost,MPI_BYTE,dest,sndtag,
				rghost,sizeof(ComGhost)*nrghost,MPI_BYTE,src,recvtag,
				MPI_COMM_WORLD,&status);
		if(myid==0 && iflag == 1){
			lpz = -nz;
			for(i=0;i<nrghost;i++) rghost[i].z += lpz;
		}
		else if(myid == nid-1 && iflag == -1){
			lpz = nz;
			for(i=0;i<nrghost;i++) rghost[i].z += lpz;
		}
		/* This is for the case when z-width < rsphere */
		if(tnneighbor > 1){
			if(myid==0) fprintf(stdout,"Now entering into too thin slab: %ld\n",tnneighbor);
			for(i=2;i<=tnneighbor;i++){
				ComGhost *addsghost;
				addsghost = NULL;
				idest = (myid+iflag*i+nid)%nid;
				isrc = (myid-iflag*i+nid)%nid;
				zwsum = 0;
				for(j=1;j<i;j++)zwsum += zw[(myid+j*iflag+nid)%nid];
				if(iflag==1){
					zmaxput = zstart+zheight;
					zminput = zstart+zheight +zwsum-rsphere;
				}
				else{
					zminput = zstart;
					zmaxput = zstart + rsphere - zwsum;
				}
				if(nneighbor >= i){
					nput = 0;
					for(j=0;j<nsghost;j++){
						if(sghost[j].z >=zminput && sghost[j].z < zmaxput)
							nput ++;
					}
					if(nput){
						addsghost = (ComGhost *)Malloc(sizeof(ComGhost)*nput,PPTR(addsghost));
						nput = 0;
						for(j=0;j<nsghost;j++){
							if(sghost[j].z >=zminput && sghost[j].z < zmaxput)
								addsghost[nput++] = sghost[j];
						}
					}
				}
				else nput =0;
				MPI_Sendrecv(&nput,1,MPI_LONG,idest,0,&nget,1,MPI_LONG,
						isrc,0,MPI_COMM_WORLD,&status);
				printf("P%d put to P%d np=%ld & get from P%d np=%ld\n"
							,myid,idest,nput,isrc,nget);
				if(nget) rghost = (ComGhost*)Realloc(rghost,sizeof(ComGhost)*(nrghost+nget));
				MPI_Sendrecv(addsghost,sizeof(ComGhost)*nput,MPI_BYTE,idest,
						10,(rghost+nrghost),sizeof(ComGhost)*nget,MPI_BYTE,isrc,
						10,MPI_COMM_WORLD,&status);
				if(myid==0 && iflag == 1){
					lpz = -nz;
					for(j=0;j<nget;j++) rghost[nrghost+j].z += lpz;
				}
				else if(myid == nid-1 && iflag == -1){
					lpz = nz;
					for(j=0;j<nget;j++) rghost[nrghost+j].z += lpz;
				}
				nrghost += nget;
				if(nput) Free(addsghost);
			}
		}

		Free(sghost);
		mbx = ceil(nx/CellWidth);
		BoundaryCell = (ComGhost **)
			Malloc(sizeof(ComGhost *)*mbx, PPTR(BoundaryCell));
		NBoundaryCell = (long *)
			Malloc(sizeof(long)*mbx, PPTR(NBoundaryCell));
		BoundaryCellLinkedList(rghost, nrghost, BoundaryCell, NBoundaryCell);
		{
			long nposmax,ngrvmax;
			nposmax = 0;
			for(i=0;i<mx;i++) nposmax = max(nposmax,NBoundaryCell[i]);
			mpos = npmaxcell * (nstride+2) *(nstride+2)*(nstride+2);
			mgrv = nposmax * nstride;
#ifdef _OPENMP
#pragma omp parallel
#endif
			{
				if(omp_get_thread_num()==0) {
					pos = (particle *)Malloc(sizeof(particle)*mpos*omp_get_num_threads(),PPTR(pos));
					grv = (TPtlStruct *) Malloc(sizeof(TPtlStruct)*mgrv*omp_get_num_threads(),PPTR(grv));
					TREECELL = (TStruct *)Malloc(sizeof(TStruct)*mgrv*omp_get_num_threads(),PPTR(TREECELL));
				}
			}
		}
#ifdef _OPENMP
#pragma omp parallel private(ix0,ix1,jxp,ngrv,npos,box,i,ix,xtran)
#endif
		{
			int idthread,nthreads;
			long myix0start,myix0end,chunk;
			particle *threadpos;
			TPtlStruct *threadgrv;
			TStruct *threadTREECELL;

			idthread = omp_get_thread_num();
			nthreads = omp_get_num_threads();
			threadpos = pos + mpos*idthread;
			threadgrv = grv + mgrv*idthread;
			threadTREECELL = TREECELL + mgrv*idthread;
			chunk = (mx+nthreads-1)/nthreads;
			myix0start = chunk*idthread;
			myix0end = min(myix0start+chunk,mx);
#ifdef _OPENMP
			for(ix0=myix0start;ix0<myix0end;ix0+=nstride){
#else
			for(ix0=0;ix0<mx;ix0+=nstride){
#endif
				xtran = -ix0*CellWidth;
				ix1 = min(ix0+nstride,mx);
				jxp = ix1-ix0+1;
				ngrv = 0;
				for(ix=ix0;ix<ix1;ix++) ngrv += NBoundaryCell[ix];
				if(ngrv > 0){
					ComGhost *tp;
					ngrv = 0;
					for(ix=ix0;ix<ix1;ix++){
						tp = BoundaryCell[ix];
						while(tp){
							threadgrv[ngrv].type = TYPE_PTL;
							threadgrv[ngrv].r[0] = tp->x + xtran;
							threadgrv[ngrv].r[1] = tp->y + ytran;
							threadgrv[ngrv].r[2] = tp->z + ztran;
#ifdef MULTIMASS
							threadgrv[ngrv].mass = tp->mass;
#else
							threadgrv[ngrv].mass = 1;
#endif
							tp = tp->next;
							ngrv ++;
						}
					}
					npos = Dump2PeriodicPosExceptZ(ix0,iy0,iz0pos,-1,jxp,-1,jyp,jzm,jzp,
							nx,ny,nz,nspace,threadpos,
							xtran,ytran,ztran);
					if(npos > 0){
#ifndef USE_GPU
						if(ngrv >= DIRECTSUM){
							box = CheckBoundingBox(ngrv,threadgrv);
							Make_Tree(threadTREECELL,threadgrv,ngrv,box,theta);
							for(i=0;i<npos;i++) treeforce(threadpos+i,theta2,threadTREECELL,threadgrv,rspheresq);
						}
						else DirectSummation(ngrv,threadgrv,npos,threadpos);
#else
						if(ngrv <= GPUCPUDIRECTSUM){
							DirectSummation(ngrv,threadgrv,npos,threadpos);
						}
						else if(ngrv <= GPUDIRECTSUM){
							Direct_Nbody(threadpos,npos, threadgrv,ngrv , GPUSPERNODE, myid,nx);
						}
						else { /* This is GPU version of tree sum */
							int ncell;
							box = CheckBoundingBox(ngrv,threadgrv);	
							ncell = Make_Tree(threadTREECELL,threadgrv,ngrv,box,theta);
							gputreeforce(threadpos,npos,threadTREECELL,ncell,threadgrv,ngrv,
								rspheresq,GPUSPERNODE,myid,nx);
						}
#endif
						PartialUpdateAccelUsingForceByParticles(vfact2,npos,threadpos,mass);
					}

				}
			} /* end of for of running index i */
		}
#ifdef _OPENMP
#pragma omp parallel
#endif
		{
				if(omp_get_thread_num()==0) {
						Free(TREECELL);
						Free(grv);
						Free(pos);
				}
		}

		Free(BoundaryCell); Free(NBoundaryCell); Free(rghost);
	}
	Free(zw);
}
