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


int nsteps[NSTEPS] ={
	167 , 191 , 213 , 235 , 256 , 276 , 295 , 314 , 332 , 350 , 368 , 418 , 434 , 450 , 466 , 481 , 496 , 511 , 526 , 541 , 556 , 570 , 584 , 598 , 612 , 626 ,
	640 , 654 , 667 , 681 , 694 , 708 , 721 , 734 , 747 , 760 , 773 , 786 , 799 , 812 , 825 , 837 , 850 , 863 , 875 , 888 , 900 , 913 , 925 , 938 , 950 , 963 , 975 , 987 , 1000,
	1012, 1024, 1037, 1049, 1061, 1074, 1086, 1098, 1110, 1123, 1135, 1147, 1160, 1172, 1184, 1197, 1209, 1221, 1234, 1246, 1258, 1271, 1283, 1296,
	1308, 1321, 1333, 1346, 1358, 1371, 1383, 1396, 1409, 1421, 1434, 1447, 1460, 1472, 1485, 1498, 1511, 1524, 1537, 1550, 1563, 1576, 1590, 1603, 1616, 1629, 1643, 1656,
	1669, 1683, 1696, 1710, 1724, 1737, 1751, 1765, 1779, 1793, 1807, 1821, 1835, 1849, 1863, 1877, 1892, 1906, 1920, 1935, 1949, 1964, 1979, 1981
};




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


IDTYPE readmergerdata(char *infile, MajorGid **smh){
	IDTYPE np,i,j,k;
	IDTYPE stnp, wtnp;
	size_t nshift = sizeof(float)*10 + sizeof(int);
	nshift = sizeof(float)*(9+NSTEPS) + sizeof(int)*NSTEPS;

	if(myid==-1){
		int kk=1;
		while(kk){
			i = i;
		}
	}

	if(myid==0){
		FILE *fp = fopen(infile,"r");
		fseek(fp, 0, SEEK_END);
		long len = (long) ftell(fp);
		fseek(fp, 0, SEEK_SET);
		fread(&boxsize,sizeof(float),1,fp);
		fread(&hubble,sizeof(float),1,fp);
		fread(&npower,sizeof(float),1,fp);
		fread(&omep,sizeof(float),1,fp);
		fread(&omepb,sizeof(float),1,fp);
		fread(&omeplam,sizeof(float),1,fp);
		fread(&bias,sizeof(float),1,fp);

		fread(&astep, sizeof(float), 1,fp);
		fread(&msteps, sizeof(int), 1,fp);
		reds = (float*)malloc(sizeof(float)*msteps);
		fread(nsteps, sizeof(int), msteps,fp);
		fread(reds, sizeof(float), msteps,fp);

		nspace = 1;

		rnx = nx;


		pmass = 3./8./M_PI/G*sqr(100.E5*hubble)*1.E6*(pc/onesolarmass);
		pmass = pmass*pow3(boxsize/hubble)/pow3(rnx/nspace)*omep;

		minmass = 30.L*pmass * 0.999L;
		printf("The minimum halo mass is %g :: %g %d %g\n",minmass, boxsize, nx, omep);

		fclose(fp);
		len = len - nshift;
		tnp = len/(sizeof(MergeHistory));
		printf("Total number of SSortMergeHistory is %ld :: %ld\n", tnp, len%(sizeof(MergeHistory)));


		for(i=0;i<msteps-1;i++){
			if(fabs(reds[i] - redshift) < 0.1 && (fabs(reds[i]-redshift) < fabs(reds[i+1]-redshift))) {
				iredshift = i;
				break;
			}
		}
		printf("The iredshift %d and redshift is %g :: %g\n",iredshift, reds[iredshift], reds[iredshift+1]);
	}
	MPI_Bcast(&boxsize, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&hubble, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&npower, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omep, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omepb, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omeplam, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&bias, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&msteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(myid!=0){
		reds = (float*)malloc(sizeof(float)*msteps);
	}
	MPI_Bcast(reds, msteps, MPI_FLOAT, 0, MPI_COMM_WORLD);

	MPI_Bcast(&iredshift, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&tnp, sizeof(IDTYPE), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&minmass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
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


	MajorGid *nbp = (MajorGid*)malloc(sizeof(MajorGid) * wtnp);
	*smh =  nbp;


	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	FILE *fp = fopen(infile,"r");
	nshift += sizeof(MergeHistory)*stnp;
	printf("P%d tries to read merger data :::: %ld\n",myid, wtnp);
	fseek(fp, nshift, SEEK_SET);

	MergeHistory *bp = (MergeHistory*)malloc(sizeof(MergeHistory)*10000);

	IDTYPE igid = gmin;
	IDTYPE icount = 0;
	IDTYPE jcount = 0;
	while((np = fread(bp,sizeof(MergeHistory),10000,fp))){
		for(i=0;i<np;i++){
			{
				nbp[icount].gid = igid;
				for(j=0;j<NSTEPS;j++) {
					if(bp[i].step[j].majorglobalid == 0 && bp[i].step[j].mbp <=0) {
						nbp[icount].majorgid[j] = -1;
					}
					else {
						nbp[icount].majorgid[j] = bp[i].step[j].majorglobalid;
					}
				}
				icount ++; 
			}
			jcount ++;

			igid ++;
			if(jcount >= wtnp ) goto out;
		}
	}
out:
	printf("P%d Total number of Icount = %ld\n",myid, icount);
	fclose(fp);
	free(bp);

	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
	return icount;
}
void Dump2Disk(FILE *fp, long shift, GidMergeHistory *bp, IDTYPE np){
	long ishift,jshift, i,j,k;
	if(np ==0) return;

	jshift = shift + bp[0].gid * sizeof(MergeHistory);
	fseek(fp, jshift, SEEK_SET);

	for(i=0;i<np;){
		IDTYPE gid = bp[i].gid;
		MergeHistory *abp = (MergeHistory*)(bp+i);

		for(j=i+1;j<np;j++) if(bp[j].gid != (bp[j-1].gid +1L)) break;
		for(k=i;k<j;k++) {
			memmove(&(abp[k-i]), &(bp[k].mergerhistory), sizeof(MergeHistory));
		}
		ishift = shift + gid*sizeof(MergeHistory);
		fseek(fp, (ishift-jshift), SEEK_CUR);
		fwrite(abp, sizeof(MergeHistory), (j-i), fp); fflush(fp);
		jshift = ishift + (j-i)*sizeof(MergeHistory);
		i  = j;
	}
}
void writehaloQ(MajorGid *sbp, IDTYPE np,char *infile, char *outfile, IDTYPE lgmin, IDTYPE lgmax){ 
	IDTYPE i,j,k; 
	int isend, iget ; 
	isend = iget =1; 
	int WGroupSize = WGROUPSIZE; 
	int src = myid-1; 
	int tgt = myid+1; 
	int itag = 18; 
	int jtag = 19; 
	MPI_Status status; 
	long shift;

	if(myid==nid-1){
		FILE *fp = fopen(outfile, "w");
		fwrite(&boxsize,sizeof(float),1,fp);
		fwrite(&hubble,sizeof(float),1,fp);
		fwrite(&npower,sizeof(float),1,fp);
		fwrite(&omep,sizeof(float),1,fp);
		fwrite(&omepb,sizeof(float),1,fp);
		fwrite(&omeplam,sizeof(float),1,fp);
		fwrite(&bias,sizeof(float),1,fp);
		fwrite(&astep, sizeof(float), 1,fp);
		fwrite(&msteps, sizeof(int), 1,fp);
		fwrite(nsteps, sizeof(int), msteps,fp);
		fwrite(reds, sizeof(float), msteps,fp);

		shift = sizeof(float)*(8+msteps) + sizeof(int)*(1+msteps);
		shift += sizeof(MergeHistory)*(lgmax-1);
		fseek(fp, shift, SEEK_SET);
		MergeHistory *bp = (MergeHistory*)malloc(sizeof(MergeHistory));
		fwrite(bp, sizeof(MergeHistory), 1, fp);
		free(bp);
		fclose(fp);

	}
	MPI_Bcast(&msteps, 1, MPI_INT, nid-1, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	shift = sizeof(float)*(8+msteps) + sizeof(int)*(1+msteps);

	MyFixedRangePqsort(sbp,np,MajorGid,IDTYPE,gid,lgmin,lgmax,MPI_COMM_WORLD);


	IDTYPE niter = 10;
	IDTYPE igrun = lgmin;
	IDTYPE nchunk = (lgmax-lgmin + niter-1)/niter;
	IDTYPE nleft = (lgmax-lgmin);
	for(i=0;i<niter;i++){

		if(myid==0){
			printf("\n\n\n ######################################################\n");
			printf(" Saving Data at iteration %ld \n",i);
			printf("\n\n\n ######################################################\n");
		}

		if(nleft< nchunk) nchunk = nleft;
		else nchunk = (lgmax-lgmin + niter-1)/niter;

		if(RANKINGROUP(myid,WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status); 
		FILE *fp = fopen(infile, "r");
		long ishift = shift + sizeof(MergeHistory)*igrun;
		GidMergeHistory *bp = (GidMergeHistory*)malloc(sizeof(GidMergeHistory)*nchunk);
		MergeHistory *abp = (MergeHistory*)bp;
		fseek(fp, ishift, SEEK_SET);
		nchunk = fread(abp, sizeof(MergeHistory), nchunk, fp);
		fclose(fp);
		if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

		for(j=nchunk-1;j>=0;j--){
			MergeHistory tmp = abp[j];
			bp[j].mergerhistory = tmp;
			bp[j].gid = igrun + j;
		}
		IDTYPE icount = 0;
		IDTYPE jj = 0;
		for(j=0;j<nchunk;j++){
restart:
			if(bp[j].gid == sbp[jj].gid){
				bp[j].gid = sbp[jj].newgid;
				for(k=0;k<NSTEPS;k++){
					bp[j].mergerhistory.step[k].majorglobalid = sbp[jj].majorgid[k];
				}
				icount ++;
				continue;
			}
			else if(bp[j].gid > sbp[jj].gid){
				jj++;
				goto restart;
			}
		}
		printf("P%d found %ld   among %ld\n",myid, icount, nchunk);
		IDTYPE ichunk = nchunk;

		MyFixedRangePqsort(bp,ichunk,GidMergeHistory,IDTYPE,gid,lgmin,lgmax,MPI_COMM_WORLD);


		if(RANKINGROUP(myid,WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,jtag,MPI_COMM_WORLD,&status); 
		FILE *wp = fopen(outfile,"r+");
		Dump2Disk(wp, shift, bp, ichunk);
		fclose(wp);
		if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) MPI_Send(&isend,1,MPI_INT,tgt,jtag,MPI_COMM_WORLD);

		free(bp);

		nleft = nleft - nchunk;
		igrun += nchunk;
	}


}
void GetLocalMinMaxofGid(MajorGid *hid,IDTYPE npnow, IDTYPE *ghmin, IDTYPE *ghmax){
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

void debug(HID *arr, IDTYPE np, int itag){
	IDTYPE i;
	for(i=0;i<np;i++){
		if(arr[i].gid == 0 ){
			printf("P%d found zero gid %ld    (%d)\n",myid,i, itag);
		}
	}
	return;
}



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


	redshift = atof(argv[1]);

	if(myid==0){
		printf("###########################################################################\n");
		printf("This program add halo quantity from FoFHaloQuantities  to the merging tree.\n");
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


	sprintf(infile,"MergingTree.dat");
	MajorGid *smbp;
	IDTYPE npsmbp = readmergerdata(infile, &smbp );

	IDTYPE gmin, gmax, lgmin, lgmax;
	GetMinMaxofGid(smbp, npsmbp, &gmin, &gmax);

	GetLocalMinMaxofGid(smbp, npsmbp, &lgmin, &lgmax);

	IDTYPE icount = 0;
	for(i=0;i<npsmbp;i++) {
		if(smbp[i].majorgid[iredshift] <0) {
			smbp[i].mgid = LONG_MAX;
			icount ++;
		}
		else smbp[i].mgid = smbp[i].majorgid[iredshift];
	}
	printf("P%d has negative halo value %ld\n",myid, icount);
	MyFixedRangePqsort(smbp,npsmbp,MajorGid,IDTYPE,mgid,lgmin,lgmax,MPI_COMM_WORLD);

	IDTYPE igrun = 0;
	int src, itag=1;
	int tgt;
	MPI_Status mpi_status;
	src = myid -1;
	tgt = myid +1;

	if(myid!=0) MPI_Recv(&igrun, 1, MPI_LONG, src, itag, MPI_COMM_WORLD, &mpi_status);
	for(i=0;i<npsmbp;i++){
		smbp[i].newgid = igrun + i;
	}
	igrun += npsmbp;
	if(myid!=nid-1) MPI_Send(&igrun, 1, MPI_LONG, tgt, itag, MPI_COMM_WORLD);

	GID *newgid = (GID*)malloc(sizeof(GID)*(npsmbp));
	IDTYPE npnewgid = npsmbp;
	for(i=0;i<npsmbp;i++){
		newgid[i].newgid = smbp[i].newgid;
		newgid[i].gid = smbp[i].gid;
	}
	MyFixedRangePqsort(newgid, npnewgid, GID, IDTYPE, gid, lgmin, lgmax, MPI_COMM_WORLD);


	jredshift = iredshift;
	/*
	goto TEST;
	*/
	for(jredshift = NSTEPS-1;jredshift>=0;jredshift--){
		IDTYPE icount = 0;
		for(i=0;i<npsmbp;i++) {
			if(smbp[i].majorgid[jredshift] <0) {
				icount ++;
				smbp[i].mgid = LONG_MAX;
			}
			else smbp[i].mgid = smbp[i].majorgid[jredshift];
		}
		printf("P%d has negative halo value %ld at ired= %d\n",myid, icount, jredshift);
		MyFixedRangePqsort(smbp,npsmbp,MajorGid,IDTYPE,mgid,lgmin,lgmax,MPI_COMM_WORLD);
		j = 0;
		for(i=0;i<npsmbp;i++){
			if(smbp[i].mgid >= LONG_MAX) continue;
recheck: 
			if(smbp[i].mgid == newgid[j].gid){
				smbp[i].majorgid[jredshift] = newgid[j].newgid;
				continue;
			}
			else if(smbp[i].mgid > newgid[j].gid){
				j++;
				if(j >= npnewgid) break;
				else goto recheck;

			}
		}
	}
TEST:
	free(newgid);
	MyFixedRangePqsort(smbp,npsmbp,MajorGid,IDTYPE,gid,lgmin,lgmax,MPI_COMM_WORLD);
	
	sprintf(outfile,"SortByMajorGidMergingTree.%.5d.dat", iredshift);
	writehaloQ(smbp, npsmbp, infile, outfile, lgmin,lgmax);


	if(myid==0) printf("Exit with success\n");
	MPI_Finalize();
	return 0;

}

