/* This code outputs the major merging tree of the most massive halos at the input time step */
/* This is version 2 which needs lesser memory space but pretty much slower. */
/* It is not complete and can not be used !!!!! 2017/06/05 */
#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<math.h>
#include<string.h>
#include<mpi.h>
#include<sys/stat.h>
#include<limits.h>

#include "merger.h"
#include "mkmerger.h"
#include "pqsort.SMG.h"
#include "SMG.h"

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

double Msun=1.989E33L;
double G= 6.672E-8L;
double pc= 3.08567802E18L;
double PI= 3.14159265354979L;
double onesolarmass=1.989E33L;

float redshift;


int nsteps[MSTEPS];




int myid, nid;
double pmass,minmass,rnx;
float boxsize,hubble,npower,omep,omepb,omeplam,bias,amax,astep,anow;
int nx,nspace;
int msteps;
float *reds;


size_t initshift;


#define MIN(a,b) ((a)<(b) ? (a): (b))
#define MAX(a,b) ((a)>(b) ? (a): (b))
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define sqr(A) ((A)*(A))
int iredshift,jredshift;
IDTYPE tnp;

void  writemergerdata(char *outfile, GidMergeHistory *bp, IDTYPE np){
	IDTYPE i,j,k;
	IDTYPE stnp, wtnp;
	int isend, iget ;
	isend = iget =1;
	int WGroupSize = WGROUPSIZE;
	int src = myid-1;
	int tgt = myid+1;
	int itag= 0;
	MPI_Status status;


	if(myid!=0)
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	FILE *wp = fopen(outfile,"a");
	printf("P%d tries to write merger data :::: %ld\n",myid, np);
	fwrite(bp, sizeof(GidMergeHistory), np, wp);
	fclose(wp);
	if(myid != nid-1) MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
}

IDTYPE readmergerdata(char *infile, GidMergeHistory **smh){
	IDTYPE np,i,j,k;
	IDTYPE stnp, wtnp;
	size_t nshift = sizeof(float)*10 + sizeof(int);
	nshift = sizeof(float)*(10+MSTEPS) + sizeof(int)*(1+MSTEPS);

	if(myid==0){
		FILE *fp = fopen(infile,"r");
		fseek(fp, 0, SEEK_END);
		long len = (long) ftell(fp);
		fseek(fp, 0, SEEK_SET);

		tnp = len/(sizeof(GidMergeHistory));
		printf("Total number of SSortMergeHistory is %ld :: %ld\n", tnp, len%(sizeof(MergeHistory)));
	}
	MPI_Bcast(&tnp, sizeof(IDTYPE), MPI_BYTE, 0, MPI_COMM_WORLD);
	wtnp = (tnp + nid-1)/nid;
	stnp = myid*wtnp;

	wtnp = MIN(wtnp, tnp - stnp);
	IDTYPE gmin = stnp;
	int isend, iget ;
	isend = iget =1;
	int WGroupSize = WGROUPSIZE;
	int src = myid-1;
	int tgt = myid+1;
	int itag= 0;
	MPI_Status status;


	GidMergeHistory *nbp = (GidMergeHistory*)malloc(sizeof(GidMergeHistory) * wtnp);
	*smh =  nbp;


	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	FILE *fp = fopen(infile,"r");
	printf("P%d tries to read merger data :::: %ld\n",myid, wtnp);
	long shift = sizeof(GidMergeHistory)*stnp;
	fseek(fp, shift, SEEK_SET);
	np = fread(nbp, sizeof(GidMergeHistory), wtnp, fp);
	fclose(fp);
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

	IDTYPE mp = 0;
	for(i=0;i<np;i++){
		if(nbp[i].gid <=0 && nbp[i].step[NSTEPS-1].mbp <=0)  {}
		else {
			nbp[mp] = nbp[i];
			mp++;
		}
	}
	np = mp;
	nbp = (GidMergeHistory*)realloc(nbp, sizeof(GidMergeHistory)*np);
	*smh = nbp;

	return np;
}
void Dump2Disk_SLOW(FILE *fp, long shift, GidMergeHistory *bp, IDTYPE np, IDTYPE lgmin, IDTYPE lgmax){
	long ishift,jshift, i,j,k;
	if(np ==0) return;

	shift = 0;


	long is, js;

	int ndiv = 80;
	IDTYPE mp = (lgmax-lgmin + ndiv - 1)/ndiv;
	GidMergeHistory *rbp = (GidMergeHistory*)malloc(sizeof(GidMergeHistory)*mp);
	fseek(fp, shift, SEEK_SET);
	is = js = 0;
	for(i=0;i<ndiv;i++){
		long ipos = ftell(fp);
		long nread = fread(rbp, sizeof(GidMergeHistory), mp, fp);
		js = is + mp;
		js = MIN(js, (lgmax-lgmin) );
		for(j=0;j<np;j++){
			if(bp[j].gid >= is && bp[j].gid<js){
				long ii = bp[j].gid -is;
				rbp[ii] = bp[j];
			}
		}
		fseek(fp, ipos, SEEK_SET);
		fwrite(rbp, sizeof(GidMergeHistory), (js-is), fp);
		is = js;
	}
	free(rbp);
}
#define Nicount (1<<16)

void AdGetLocalMinMaxofGid(GidMergeHistory *hid,IDTYPE npnow, IDTYPE *ghmin, IDTYPE *ghmax, IDTYPE *Gmin, IDTYPE *Gmax){
	IDTYPE i,j,k;
	IDTYPE rmin,rmax,tmin,tmax;
	rmin = INT_MAX;
	rmax = INT_MIN;
	IDTYPE icount[Nicount];
	IDTYPE jcount[Nicount];
	IDTYPE gmin,gmax;

	for(i=0;i<npnow;i++){
		rmin = MIN(rmin, hid[i].gid);
		rmax = MAX(rmax, hid[i].gid);
	}
	MPI_Allreduce(&rmin, &gmin, 1, MPI_LONG, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&rmax, &gmax, 1, MPI_LONG, MPI_MAX, MPI_COMM_WORLD);
	printf("P%d has local min/max %ld  %ld among %ld %ld \n", myid, rmin, rmax, gmin, gmax);
	*Gmin = gmin;
	gmax ++;
	*Gmax = gmax;




	IDTYPE ispacing = (gmax-gmin+Nicount)/Nicount;
	IDTYPE ioffset=0,joffset,nTot, lTot;
	MPI_Status mpi_status;
	if(myid!= 0) MPI_Recv(&ioffset, 1, MPI_LONG, myid-1, 0, MPI_COMM_WORLD, &mpi_status);
	joffset = ioffset + npnow;
	if(myid!=nid-1) MPI_Send(&joffset, 1, MPI_LONG, myid+1, 0, MPI_COMM_WORLD);

	MPI_Reduce(&npnow, &nTot, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&nTot, 1, MPI_LONG,  0, MPI_COMM_WORLD);

	for(i=0;i<Nicount;i++) icount[i] = 0;
	for(i=0;i<npnow;i++) {
		int j;
		if(hid[i].gid == LONG_MAX) { 
			j = (double)(i + ioffset)/(double)nTot * Nicount;
		}
		else {
			j = (hid[i].gid-gmin)/ispacing;
		}
		if(j<Nicount && j >=0){
			icount[j] ++;
		}
	}
	MPI_Reduce(icount, jcount, Nicount, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(jcount, Nicount, MPI_LONG, 0, MPI_COMM_WORLD);
	lTot = (nTot + nid-1)/nid;
	IDTYPE isum = 0;
	IDTYPE irun = 0;
	if(myid==0) *ghmin = gmin;
	for(i=0;i<Nicount;i++){
		isum += jcount[i];
		if(isum >= lTot) {
			isum = isum - lTot;
			irun ++ ;
			if(irun == myid+1){
				*ghmax = (i+1)*ispacing + gmin;
			}
			else if(irun == myid){
				*ghmin = (i+1)*ispacing + gmin;
			}
		}
	}
	if(myid== nid-1) *ghmax = gmax + 1;

	printf("P%d local data min/max = %ld %ld\n",myid, *ghmin, *ghmax);
	/*
	IDTYPE gstep;
	IDTYPE istep = (gmax-gmin + nid - 1)/nid;
	IDTYPE jrun= 0;
	for(i=0;i<nid;i++){
		if(i < nid/2) {
			gstep = istep * (0.3L + 0.7L*(i/(nid/2.L)));
		}
		else {
			gstep = (gmax-gmin-jrun)/(nid-i);
		}
		if(myid==i){
			*ghmin = jrun;
			*ghmax = jrun+gstep;
			if(myid == nid-1) *ghmax = gmax+1;
			return;
		}
		jrun += gstep;
	}
	return;
	*/
}
void Dump2Disk_OLD(FILE *fp, long shift, GidMergeHistory *bp, IDTYPE np, IDTYPE lgmin, IDTYPE lgmax){
	long ishift,jshift, i,j,k;
	if(np ==0) return;
	shift = 0;
	jshift = shift + (bp[0].gid-lgmin) * sizeof(GidMergeHistory);
	fseek(fp, jshift, SEEK_SET);
	for(i=0;i<np;){
		IDTYPE gid = bp[i].gid;
		for(j=i+1;j<np;j++) if(bp[j].gid != (bp[j-1].gid +1L)) break;
		ishift = shift + (gid-lgmin)*sizeof(GidMergeHistory);
		fseek(fp, (ishift-jshift), SEEK_CUR);
		fwrite((bp+i), sizeof(GidMergeHistory), (j-i), fp); 
		jshift = ishift + (j-i)*sizeof(GidMergeHistory);
		i  = j;
	}
}
void Dump2Disk(FILE *fp, long shift, GidMergeHistory *bp, IDTYPE np, IDTYPE lgmin, IDTYPE lgmax){
	long ishift,jshift, i,j,k;
	if(np ==0) return;
	fseek(fp, 0L, SEEK_END);
	fwrite(bp, sizeof(GidMergeHistory), np,fp);
}

void GetLocalMinMaxofGid(GidMergeHistory *hid,IDTYPE npnow, IDTYPE *ghmin, IDTYPE *ghmax){
	IDTYPE i,j,k;
	IDTYPE rmin,rmax,tmin,tmax;
	rmin = INT_MAX;
	rmax = INT_MIN;
	for(i=0;i<npnow;i++){
		if(rmin > hid[i].gid) rmin = hid[i].gid;
		if(rmax < hid[i].gid) rmax = hid[i].gid;
	}
	printf("P%d local data min/max = %ld %ld\n",myid, rmin,rmax);
	*ghmin = rmin;
	*ghmax = rmax+1;
	return;
}

void GetMinMaxofGid(MajorGid *hid,IDTYPE npnow, IDTYPE *gmin, IDTYPE *gmax){
	IDTYPE i,j,k;
	IDTYPE rmin,rmax,tmin,tmax;
	rmin = INT_MAX;
	rmax = INT_MIN;
	for(i=0;i<npnow;i++){
		if(rmin > hid[i].gid) rmin = hid[i].gid;
		if(rmax < hid[i].gid) rmax = hid[i].gid;
	}
	printf("P%d local data min/max = %ld %ld\n",myid, rmin,rmax);
	long lrmin,ltmin,lrmax,ltmax;
	lrmin = rmin;
	lrmax = rmax;
	MPI_Reduce(&lrmin,&ltmin, 1, MPI_LONG,MPI_MIN,0, MPI_COMM_WORLD);
	MPI_Reduce(&lrmax,&ltmax, 1, MPI_LONG,MPI_MAX,0, MPI_COMM_WORLD);
	if(myid==0) { rmin = ltmin; rmax=ltmax;}
	MPI_Bcast(&rmin, 1, MPI_IDTYPE, 0,MPI_COMM_WORLD);
	MPI_Bcast(&rmax, 1, MPI_IDTYPE, 0,MPI_COMM_WORLD);
	*gmin = rmin;
	*gmax = rmax;
}

/*
void debug(HID *arr, IDTYPE np, int itag){
	IDTYPE i;
	for(i=0;i<np;i++){
		if(arr[i].gid == 0 ){
			printf("P%d found zero gid %ld    (%d)\n",myid,i, itag);
		}
	}
	return;
}
*/



int main(int argc, char **argv){
	IDTYPE i,j,k;
	IDTYPE np1,np2,ip;
	TrHalo *nmbp,imbp;
	MbpRV imbprv;
	FILE *fp1;
	HidMass *hidmass;
	char infile[190],outfile[190];


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);


	int nfiles = atof(argv[2]);

	if(myid==0){
		printf("###########################################################################\n");
		printf("This program writes down the sorted mbp data.                              \n");
		printf("###########################################################################\n");
	} 
	{ 
		char hostname[190]; 
		hostname[189] = '\0'; 
		gethostname(hostname,189); 
		printf("P%d host name is %s with pid= %d\n",myid,hostname, getpid()); 
	}


	FILE *fp;
	float nredshift;
	/*
	IDTYPE npmbp;
	SortMergeHistory *mbp;
	*/
	MPI_Status status;


	if(myid==-1){
		int kk=1;
		while(kk){
			i = i;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);

	sprintf(outfile,"%s", argv[1]);

	for(i=0;i<nfiles;i++){
		IDTYPE lgmin, lgmax;
		IDTYPE gmin, gmax;
		sprintf(infile,"%s.%.5d",argv[1],i);
		GidMergeHistory *smbp;
		IDTYPE npsmbp = readmergerdata(infile, &smbp);

		AdGetLocalMinMaxofGid(smbp, npsmbp, &lgmin, &lgmax, &gmin, &gmax);

		MyFixedRangePqsort(smbp,npsmbp,GidMergeHistory,IDTYPE,gid,lgmin,lgmax,MPI_COMM_WORLD);

		writemergerdata(outfile, smbp, npsmbp);
		free(smbp);

	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myid==0) printf("Exit with success\n");
	MPI_Finalize();
	return 0;

}

