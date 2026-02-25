/* pgcc -o repartition repartition.c -lm -g header.c -DINDEX */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#define DEFINE_SIM_PARA
#	include "pmheader.h"
#undef DEFINE_SIM_PARA
#include "params.h"

int jwrite(char *, SimParameters ,int ,pmparticletype *, long );
SimParameters jread(char *, long , long *);
void splitmigrate(pmparticletype *,long ,float );

float zsplit;

int sorttype(const void *a,const void *b){
	pmparticletype *aa,*bb;
	aa = (pmparticletype *)a;
	bb = (pmparticletype *)b;
	if(ZofP(aa) < ZofP(bb)) return -1;
	else if(ZofP(aa) > ZofP(bb)) return +1;
	else return 0;
}
pmparticletype *pmparticles;

int main(int argc, char **argv){
	long np,nfile,mfile,mp,ifile,iread;
	long i,j,k;
	char insyncfile[190],outsyncheader[190];
	char insyncheader[190],parafile[190];
	FILE *fp,*wp;
	SimParameters outpar,readpar;
	nfile = atoi(argv[2]);
	sprintf(outsyncheader,"%s",argv[3]);
	pmparticles = (pmparticletype*)malloc(sizeof(pmparticletype)*100);
	{
		sprintf(parafile,"%s00000",argv[1]);
		fp = fopen(parafile,"r");
		readpar=read_head(fp);
		mfile = atoi(argv[4]);
		mp = (long)readpar.nx*(long)readpar.ny*(long)readpar.nz 
				/(readpar.nspace*readpar.nspace*readpar.nspace);
		mp = mp / mfile;
		fclose(fp);
	}


	np = 0;
	iread = ifile = 0;
	for(i=0;i<mfile;i++){
		while(np < mp){
			long oldnp;
			sprintf(insyncfile,"%s%.5d",argv[1],iread);
			oldnp = np;
			readpar = jread(insyncfile,np, &np);
			qsort(pmparticles+oldnp,(np-oldnp),sizeof(pmparticletype),sorttype);
			iread ++;
		}	
		outpar = readpar;

		if(ifile == 0) outpar.zmin = 0;
		else outpar.zmin = ZofP(pmparticles);

		if(ifile == mfile-1) outpar.zmax = simpar.nz;
		else {
			if(np > mp) outpar.zmax = ZofP(pmparticles+mp);
			else if(np ==mp) outpar.zmax = readpar.zmax;
		}

		if(ifile == mfile-1) outpar.np = np;
		else outpar.np = mp;

		outpar.nid = mfile;
		outpar.myid = ifile;

		jwrite(outsyncheader, outpar,ifile,pmparticles,outpar.np);
		ifile ++;
		for(j=0;j<np-outpar.np;j++){
			pmparticles[j] = pmparticles[j+outpar.np];
		}
		np = np - outpar.np;
	}

}
SimParameters  jread(char *rvfilename, long mp, long *newnp)
{
	FILE *rvfile;
	int i, j, k, count;
	long nmax,np;
	int nsend, nrecv;
	long NTotal, NLocal, NLocal0, NTotalcheck;
        int Ntotal;
	long nbytes;
	float m, mrv[6];
	char rvname[80],filename[80];
	float ksum, totalken;
	pmparticletype *p;
	int src,tgt,WGroupSize;
	int isend,iget,itag=1;
	SimParameters readpar;
	sprintf(filename,"%s%",rvfilename);
	{
		rvfile = fopen(filename,"r");
		readpar=read_head(rvfile);
		readpar.local_z_start = rint(readpar.zmin);
		np = readpar.np;
		pmparticles = (pmparticletype*)realloc(pmparticles,sizeof(pmparticletype)*(mp+np));
       	fprintf(stderr,"starts reading data. %ld",np);
		fflush(stderr);
		fread(pmparticles+mp,sizeof(pmparticletype),np,rvfile);
       	fprintf(stderr," --  finishes reading data, %ld \n", np);
		fflush(stderr);
		*newnp = mp+np;
	}
	return readpar;
}
int jwrite(char *rvprefix, SimParameters outpar,int myid,pmparticletype *localpmparticles,
				long mp) {
	FILE *fp;
	int i,j,k, p;
	int count;
	int ncount, nblock, nremain, ncheck;
	int dummy=0;
	float m;
	char *fileheader = FILEHEADER;
	char procname[3],number[3];
	char filename[80];
	size_t nwritecheck;
	int oldstepcount;
	int src,tgt,WGroupSize;
	int isend,iget,itag=1;
	sprintf(outpar.rvfilename,"%s",rvprefix);

	{
		long wmp;
		sprintf(filename,"%s%.5d",rvprefix,myid);
		fprintf(stderr," writing file: %s %d ",filename,mp);
		fp = fopen(filename,"w");
		write_head(fp,outpar);
		ncheck = mp;
		count = 0;
		wmp = fwrite(localpmparticles,sizeof(pmparticletype),mp,fp);
		fclose(fp);
		fprintf(stderr," written: %ld\n", wmp);
	}
}
void splitmigrate(pmparticletype *pmparticles,long np,float zsplit){
  	pmparticletype *p,tmp,*ptr,*bp,*sp;
	pmparticletype *sndp,*lp,*rp,*bp2send;
	long nrecv,nsend,commsize;
	long i,j,k,sndnp;
	float zmax;
#ifdef XYZDBL
	double zp,zup,zdown;
#else
	float zp,zup,zdown;
#endif
	long niter,tniter;
	sp = pmparticles;
	{
		rp = pmparticles + np;
		for(lp=sp;lp<rp;){
			zp = ZofP(lp);
			if(zp >= zsplit){
				rp --;
				tmp = *lp;
				*lp = *rp;
				*rp = tmp;
			}
			else lp++;
		}

	}
}
