/* This is the first phase of the parallel FoF finding algorithm */
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
#include "params.h"

#define DEFINE_PRE_FOF
#include "prefof.h"
#undef DEFINE_PRE_FOF
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
float maintreetime;
long mx,my,mz;
treeparticletype **BasicCell;
long *NBasicCell;
#ifndef _OPENMP
#define omp_get_thread_num() nullfct0()
#define omp_get_num_threads() nullfct1()
#endif


int NTotal;


int dest,src;
int idest,isrc;
int myid,nid;
void FoFParticleCellLinkedList(treeparticletype *particles, long np, float zstart,
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
		ix = min(ix,mxm1); iy = min(iy,mym1); iz = min(iz,mzm1);
		mpixel = ix+mx*(iy+my*iz);
		NBasicCell[mpixel] ++;
		p = BasicCell[mpixel];
		BasicCell[mpixel] = tp;
		tp->next = p;
		tp++;
	}
}
FoFBox CheckBoundingFoFBox(long npoint,FoFTPtlStruct *point){
	float xmin,xmax,ymin,ymax,zmin,zmax;
	float width;
	long i;
	int si;
	FoFBox box;
	xmin = ymin = zmin = 1.e25;
	xmax = ymax = zmax = -1.e25;
	for(i=0;i<npoint;i++){
		xmin = min(xmin,point[i].r[0]);
		ymin = min(ymin,point[i].r[1]);
		zmin = min(zmin,point[i].r[2]);
		xmax = max(xmax,point[i].r[0]);
		ymax = max(ymax,point[i].r[1]);
		zmax = max(zmax,point[i].r[2]);
	}
	width = xmax-xmin;
	width = max(width,ymax-ymin);
	width = max(width,zmax-zmin);
	if(width !=0){
		box.xi = xmin-width*1.e-6;
		box.yi = ymin-width*1.e-6;
		box.zi = zmin-width*1.e-6;
		box.xf = xmax+width*1.e-6;
		box.yf = ymax+width*1.e-6;
		box.zf = zmax+width*1.e-6;
		box.maxwidth = width*(1.+2.e-6);
	}
	else {
		box.xi = xmin-1.e-5;
		box.yi = ymin-1.e-5;
		box.zi = zmin-1.e-5;
		box.xf = xmax+1.e-5;
		box.yf = ymax+1.e-5;
		box.zf = zmax+1.e-5;
		box.maxwidth = 2.e-5;
	}
	return box;
}
long Dump2FoFPtl(long ix0,long iy0,long iz0,long ixm,long ixp,
		long iym,long iyp,long izm,long izp,FoFTPtlStruct *point, 
		int nx,int ny,int nz,int nspace, double xtran, double ytran,
		double ztran){
	long ix,iy,iz,mcell;
	long ixi,iyi,izi,ixf,iyf,izf;
	long ntmp,npoint;
	long zshift,yzshift,mxy;
	treeparticletype *p;
	npoint = 0;
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
					point[npoint].type = TYPE_PTL;
					point[npoint].r[0] = XofP(p)+xtran;
					point[npoint].r[1] = YofP(p)+ytran;
					point[npoint].r[2] = ZofP(p)+ztran;
					point[npoint].bp = p;
					/*
					point[npoint].sibling = &(point[npoint+1]);
					*/
					point[npoint].haloindx = -1;
					p = p->next;
					npoint++;
				}
			}
		}
	}
	/*
	if(npoint>0) point[npoint-1].sibling = NULL;
	*/
	return npoint;
}
int CountHaloCandidates(treeparticletype *particles,int nhalo,HaloBound *halo,
		FoFTPtlStruct *threadpoint,int npoint, float  fof_link, FoFBox box){
	int i,mh;
	FoFTPtlStruct *tmp;
	int localnp;
	for(i=0;i<nhalo;i++){
		halo[i].xmin = 1.E20;
		halo[i].xmax =-1.E20;
		halo[i].ymin = 1.E20;
		halo[i].ymax =-1.E20;
		halo[i].zmin = 1.E20;
		halo[i].zmax =-1.E20;
		halo[i].dumpflag = 0;
		halo[i].sibling = NULL;
		halo[i].nmem = 0;
	}
	for(i=0;i<npoint;i++){
		mh = threadpoint[i].haloindx;
		halo[mh].nmem ++;
		tmp = halo[mh].sibling;
		halo[mh].sibling = &(threadpoint[i]);
		threadpoint[i].sibling = tmp;
		halo[mh].xmin = min(halo[mh].xmin,threadpoint[i].r[0]);
		halo[mh].xmax = max(halo[mh].xmax,threadpoint[i].r[0]);
		halo[mh].ymin = min(halo[mh].ymin,threadpoint[i].r[1]);
		halo[mh].ymax = max(halo[mh].ymax,threadpoint[i].r[1]);
		halo[mh].zmin = min(halo[mh].zmin,threadpoint[i].r[2]);
		halo[mh].zmax = max(halo[mh].zmax,threadpoint[i].r[2]);
	}
	localnp = 0;
	for(i=0;i<nhalo;i++){
		if(halo[i].xmin > fof_link && halo[i].ymin > fof_link && 
				halo[i].zmin > fof_link && halo[i].xmax < box.xf - fof_link && 
				halo[i].ymax < box.yf - fof_link && halo[i].zmax < box.zf - fof_link) {
			if(halo[i].nmem >= MINFOFNP) {
				halo[i].dumpflag = 1;
				localnp += halo[i].nmem;
			}
		}
		else{
			halo[i].dumpflag = 1;
			localnp += halo[i].nmem;
		}
	}
	for(i=0;i<nhalo;i++){
		if(halo[i].dumpflag){
			tmp = halo[i].sibling;
			while(tmp){
				int mp;
				mp = (tmp->bp)-particles;
				SET_FoFflag((mp/NBitsperByte),(mp%NBitsperByte));
				tmp = tmp->sibling;
			}
		}
	}
	return localnp;
}
#define NPSAVE (1<<20)
#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))
void Dump2Disk(treeparticletype *particles, long np){
	long i,j,k,nsave;
	int jid;
	pmparticletype *pbuff;
	FILE *wp;
	char outfile[190];
	int myid,src,tgt,nid;
	int WGroupSize,isend,iget,itag=100;
	MPI_Status status;
	SimParameters nowsimpar;

	long npsave,nporg;
	myid = simpar.myid;
	nid = simpar.nid;

	nowsimpar = simpar;
	npsave = 0;
	for(i=0;i<np;i++){
		if(IS_FoFflag((i/NBitsperByte),(i%NBitsperByte))) npsave++;
	}
	nowsimpar.np = npsave;

	WGroupSize = WGROUPSIZE;
	src = myid-1;
	tgt = myid+1;
	if(RANKINGROUP(myid,WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);


	sprintf(outfile,"PreFoF.%.5d%.5d",nowsimpar.stepcount,myid);
	wp = fopen(outfile,"w");
	write_head(wp,nowsimpar);


	pbuff = (pmparticletype *)Malloc(sizeof(pmparticletype)*NPSAVE,PPTR(pbuff));
	nsave = 0;
	for(i=0;i<np;i++){
		if(IS_FoFflag((i/NBitsperByte),(i%NBitsperByte))){
			pbuff[nsave].x = particles[i].x;
			pbuff[nsave].y = particles[i].y;
			pbuff[nsave].z = particles[i].z;
			pbuff[nsave].vx = particles[i].vx;
			pbuff[nsave].vy = particles[i].vy;
			pbuff[nsave].vz = particles[i].vz;
			pbuff[nsave].indx = particles[i].indx;
			nsave ++;
			if(nsave == NPSAVE){
				fwrite(pbuff,sizeof(pmparticletype),nsave,wp);
				nsave = 0;
			}
		}
	}

	if(nsave>0) fwrite(pbuff,sizeof(pmparticletype),nsave,wp);
	fclose(wp);
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
	Free(pbuff);
}
#undef RANKINGROUP
#undef GROUPID
#undef NPSAVE

float PreFoF(treeparticletype *particles,long np){
	float mass,zstart,zfinal,zheight;
	int nx,ny,nz,nspace;
	long i,j,k;
	long ix0,iy0,iz0;
	long ix2,iy2,iz2;
	long jxp,jyp,jzp;
	long ncell,mcell,mp;
	long npos;
	long mpos,maxnpincell;
	HaloBound *halo;
	particle *pos;
	FoFTPtlStruct *point;
	FoFTStruct *TREECELL;
	treeparticletype *p,*pdump;
	FoFBox box;
	double xtran,ytran,ztran;
	float fof_link=0.2;
	particle *linked;
	INT8  nmaxdump,*localnpdump;
	long npdump=0;
#ifdef _OPENMP
	int si; /* for integer*4 value of loop index in OpenMP*/
#endif
	mass = simpar.nspace*simpar.nspace*simpar.nspace;
	zstart = simpar.zmin;
	zfinal = simpar.zmax;
	zheight = (simpar.zmax-simpar.zmin);
	nx = simpar.nx;
	ny = simpar.ny;
	nz = simpar.nz;
	nspace = simpar.nspace;
	myid = simpar.myid;
	nid = simpar.nid;

	nfofflag = (np+NBitsperByte-1)/NBitsperByte;
	fofflag = (unsigned char *)Malloc(nfofflag,PPTR(fofflag));
	for(i=0;i<nfofflag;i++) fofflag[i] = FoF_RESET;

	NTotal = 0;



	mx = (long)ceil((double)nx/(double)CellWidth);
	my = (long)ceil((double)ny/(double)CellWidth);
	mz = (long)ceil((double)zheight/(double)CellWidth);

	p = particles;
	BasicCell = (treeparticletype **) 
		Malloc(sizeof(treeparticletype *)*mx*my*mz,PPTR(BasicCell));
	NBasicCell = (long *) Malloc(sizeof(long)*mx*my*mz,PPTR(NBasicCell));
	FoFParticleCellLinkedList(particles,np,zstart,nx,ny,nz,nspace);
	{
		long nposmax,npointmax;
		nposmax = 0;
		for(i=0;i<mx*my*mz;i++) nposmax = max(nposmax,NBasicCell[i]);
		maxnpincell = nposmax * nstride*nstride*nstride;

#ifdef _OPENMP
#pragma omp parallel
#endif
		{
			if(omp_get_thread_num()==0) {
				point = (FoFTPtlStruct *) Malloc(sizeof(FoFTPtlStruct)*maxnpincell*
						omp_get_num_threads(),PPTR(point));
				TREECELL = (FoFTStruct *)Malloc(sizeof(FoFTStruct)*maxnpincell*
						omp_get_num_threads(),PPTR(TREECELL));
				linked = (particle *)Malloc(sizeof(particle)*maxnpincell*
						omp_get_num_threads(),PPTR(linked));
				halo = (HaloBound*)Malloc(sizeof(HaloBound)*maxnpincell*
						omp_get_num_threads(),PPTR(halo));
				localnpdump = (INT8 *)Malloc(sizeof(INT8)*maxnpincell*
						omp_get_num_threads(),PPTR(localnpdump));
				/*
				nmaxdump = freespace()/sizeof(treeparticletype)-1000L;
				pdump = (treeparticletype*)Malloc(sizeof(treeparticletype)*
						nmaxdump,PPTR(pdump));
						*/
			}
		}
	}
	{
		box.xi = box.yi = box.zi = 0;
		box.xf = box.yf = box.zf = CellWidth*nstride;
		box.maxwidth = CellWidth*nstride;
	}


	for(iz0=0;iz0<mz;iz0+=nstride){
		jzp = min(iz0+nstride,mz)-iz0+1;
		/* This OpenMP struct could have the data racing conditions 
		 * if "chunk" = 1 and 2 or possibly 3 (?).
		 * And "my" should be divisible by "nthreads" and */
#ifdef _OPENMP
#pragma omp parallel private(iy0,jyp,jxp,ix0,i,box,xtran,ytran,ztran)
#endif
		{
			long npoint;
			int siy0;
			int threadID,nthreads;
			long myiy0start,myiy0end,chunk;
			FoFTPtlStruct *threadpoint;
			FoFTStruct *threadTREECELL;
			particle *threadlinked;
			HaloBound *threadhalo;
			threadID = omp_get_thread_num();
			nthreads = omp_get_num_threads();
			threadpoint = point + maxnpincell*threadID;
			threadTREECELL = TREECELL + maxnpincell*threadID;
			threadlinked = linked + maxnpincell*threadID;
			threadhalo = halo + maxnpincell*threadID;

			box.xi = box.yi = box.zi = 0;
			ztran = -(iz0*CellWidth+zstart);
			box.zf = min(CellWidth*nstride,(zfinal+ztran));
			box.maxwidth = CellWidth*nstride;

			chunk = (my+nthreads-1)/nthreads;
			myiy0start = chunk*threadID;
			myiy0end = min(myiy0start+chunk,my);
#ifdef _OPENMP
			for(iy0=myiy0start;iy0<myiy0end;iy0+=nstride)
#else
			for(iy0=0;iy0<my;iy0+=nstride)
#endif
			{
				ytran = -iy0*CellWidth;
				box.yf = min(CellWidth*nstride,(ny+ytran));
				jyp = min(iy0+nstride,my)-iy0+1;
				for(ix0=0;ix0<mx;ix0+=nstride){
					int nhalo=0;
					xtran = -ix0*CellWidth;
					box.xf = min(CellWidth*nstride,(nx+xtran));
					jxp = min(ix0+nstride,mx)-ix0+1;
					npoint = Dump2FoFPtl(ix0,iy0,iz0,0,nstride,0,nstride,0,nstride,
							threadpoint,nx,ny,nz,nspace,xtran,ytran,ztran);
					if(npoint>0){
						/*
						printf("P%d has box %g %g %g :%g\n",myid,box.xf,box.yf,
										box.zf,box.maxwidth);
						*/
						FoF_Make_Tree(threadTREECELL,threadpoint,npoint,box);
						for(i=0;i<npoint;i++){
							if(threadpoint[i].included == NO){
								int num;
								particle p;
								p.x = threadpoint[i].r[0];
								p.y = threadpoint[i].r[1];
								p.z = threadpoint[i].r[2];
								num=pnew_fof_link(&p,fof_link,threadTREECELL,
										threadpoint,threadlinked,nhalo);
								nhalo ++;
							}
						}
						localnpdump[threadID] = 
							CountHaloCandidates(particles,nhalo,threadhalo,threadpoint,npoint,
									fof_link, box);
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
				/*
				Free(pdump);
				*/
				Free(localnpdump);
				Free(halo);
				Free(linked);
				Free(TREECELL);
				Free(point);
			}
	}
	Free(NBasicCell);
	Free(BasicCell);

	Dump2Disk(particles,np);

	Free(fofflag);

	return maintreetime;
}

