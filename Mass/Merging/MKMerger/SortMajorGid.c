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

/*
int nsteps[NSTEPS] ={
	167 , 191 , 213 , 235 , 256 , 276 , 295 , 314 , 332 , 350 , 368 , 418 , 434 , 450 , 466 , 481 , 496 , 511 , 526 , 541 , 556 , 570 , 584 , 598 , 612 , 626 ,
	640 , 654 , 667 , 681 , 694 , 708 , 721 , 734 , 747 , 760 , 773 , 786 , 799 , 812 , 825 , 837 , 850 , 863 , 875 , 888 , 900 , 913 , 925 , 938 , 950 , 963 , 975 , 987 , 1000,
	1012, 1024, 1037, 1049, 1061, 1074, 1086, 1098, 1110, 1123, 1135, 1147, 1160, 1172, 1184, 1197, 1209, 1221, 1234, 1246, 1258, 1271, 1283, 1296,
	1308, 1321, 1333, 1346, 1358, 1371, 1383, 1396, 1409, 1421, 1434, 1447, 1460, 1472, 1485, 1498, 1511, 1524, 1537, 1550, 1563, 1576, 1590, 1603, 1616, 1629, 1643, 1656,
	1669, 1683, 1696, 1710, 1724, 1737, 1751, 1765, 1779, 1793, 1807, 1821, 1835, 1849, 1863, 1877, 1892, 1906, 1920, 1935, 1949, 1964, 1979, 1981
};
*/




int myid, nid;
double pmass,minmass,rnx;
float boxsize,hubble,npower,omep,omepb,omeplam,bias,amax,astep,anow, wlam0, wlam1;
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
	nshift = sizeof(float)*(10+MSTEPS) + sizeof(int)*(1+MSTEPS);

	if(myid==-1){
		int kk=1;
		while(kk){
			i = i;
		}
	}

	if(myid==0){
		printf("P%d is now opening %s\n", myid, infile);fflush(stdout);
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
		fread(&wlam0,sizeof(float),1,fp);
		fread(&wlam1,sizeof(float),1,fp);
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
		iredshift = NSTEPS-1;
		printf("The iredshift %d and redshift is %g :: %g\n",iredshift, reds[iredshift], reds[iredshift+1]);
	}

	/* MSTEPS should be input integer variable */
	msteps = MSTEPS;

	MPI_Bcast(&boxsize, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&hubble, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&npower, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omep, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omepb, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&omeplam, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&wlam0, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&wlam1, 1, MPI_FLOAT, 0, MPI_COMM_WORLD);
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
	char outfile1[100];
	sprintf(outfile1,"%s.%.5d", outfile, myid);
	{
		FILE *wp = fopen(outfile1,"w");
		fclose(wp);
	}

	if(myid==0){
		FILE *fp = fopen(outfile, "w");
		fwrite(&boxsize,sizeof(float),1,fp);
		fwrite(&hubble,sizeof(float),1,fp);
		fwrite(&npower,sizeof(float),1,fp);
		fwrite(&omep,sizeof(float),1,fp);
		fwrite(&omepb,sizeof(float),1,fp);
		fwrite(&omeplam,sizeof(float),1,fp);
		fwrite(&wlam0,sizeof(float),1,fp);
		fwrite(&wlam1,sizeof(float),1,fp);
		fwrite(&bias,sizeof(float),1,fp);
		fwrite(&astep, sizeof(float), 1,fp);
		fwrite(&msteps, sizeof(int), 1,fp);
		fwrite(nsteps, sizeof(int), msteps,fp);
		fwrite(reds, sizeof(float), msteps,fp);
		fclose(fp);
	}
	MPI_Bcast(&msteps, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	shift = sizeof(float)*(10+msteps) + sizeof(int)*(1+msteps);

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
		FILE *fpsigv = fopen("MergerSgiv.dat", "r");
		long ishift = shift + sizeof(MergeHistory)*igrun;
		long isigvshift = sizeof(SigV)*igrun;
		MergeHistory *abp = (MergeHistory*)malloc(sizeof(GidMergeHistory)*nchunk);
		SigV *sigv = (SigV*)malloc(sizeof(SigV)*nchunk);
		GidMergeHistory *bp = (GidMergeHistory*)abp;

		fseek(fp, ishift, SEEK_SET);
		fseek(fpsigv, isigvshift, SEEK_SET);
		nchunk = fread(abp, sizeof(MergeHistory), nchunk, fp);
		long mchunk = fread(sigv, sizeof(SigV), nchunk, fpsigv);
		fclose(fp);
		fclose(fpsigv);
		printf("P%d passed reading data   %ld : %ld\n",myid, nchunk, mchunk);
		if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

		for(j=nchunk-1;j>=0;j--){
			MergeHistory tmp = abp[j];
			for(k=0;k<NSTEPS;k++) bp[j].step[k] = tmp.step[k];
			bp[j].sigv = sigv[j];
			bp[j].gid = igrun + j;
		}
		printf("P%d passed inserting data\n",myid);

		free(sigv);

		IDTYPE jj=0;
		for(j=0;j<nchunk;j++){
			int ll = 0;
			for(k=0;k<NSTEPS;k++) {
				if(bp[j].step[k].majorglobalid <= 0 && bp[j].step[k].mbp <=0) {
					ll ++ ;
				}
			}
			if(ll < NSTEPS){
				bp[jj] = bp[j];
				jj++;
			}
		}
		nchunk = jj;

		bp = (GidMergeHistory*)realloc(bp, sizeof(GidMergeHistory)*nchunk);
		IDTYPE icount = 0;
		jj = 0;
		for(j=0;j<nchunk;j++){
restart:
			if(bp[j].gid == sbp[jj].gid){
				bp[j].gid = sbp[jj].newgid;
				for(k=0;k<NSTEPS;k++){
					bp[j].step[k].majorglobalid = sbp[jj].majorgid[k];
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

		printf("P%d after Sorting\n",myid);

		if(RANKINGROUP(myid,WGroupSize) != 0 ) MPI_Recv(&iget,1,MPI_INT,src,jtag,MPI_COMM_WORLD,&status); 
		FILE *wp = fopen(outfile1,"r+");
		Dump2Disk(wp, shift, bp, ichunk, lgmin, lgmax);
		fclose(wp);
		if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) MPI_Send(&isend,1,MPI_INT,tgt,jtag,MPI_COMM_WORLD);

		free(bp);

		nleft = nleft - nchunk;
		igrun += nchunk;
	}


}
#define Nicount (1<<16)
void HAMB(int myid){ 
	char hostname[190]; 
	gethostname(hostname,190); 
	printf("P%d is on %s with pid=%d\n",myid,hostname,getpid()); 
	if(myid == 0){ 
		long kkk = 1; 
		while(kkk) { 
			kkk = 1; 
		} 
	}
}


void AdGetLocalMinMaxofGid(MajorGid *hid,IDTYPE npnow, IDTYPE *ghmin, IDTYPE *ghmax, IDTYPE gmin, IDTYPE gmax){
	IDTYPE i,j,k;
	IDTYPE rmin,rmax,tmin,tmax;
	rmin = INT_MAX;
	rmax = INT_MIN;
	IDTYPE icount[Nicount];
	IDTYPE jcount[Nicount];
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
		if(hid[i].mgid == LONG_MAX) { 
			j = (double)(i + ioffset)/(double)nTot * Nicount;
		}
		else {
			j = (hid[i].mgid-gmin)/ispacing;
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
void GetLocalMinMaxofGid(MajorGid *hid,IDTYPE npnow, IDTYPE *ghmin, IDTYPE *ghmax, IDTYPE gmin, IDTYPE gmax){
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


	redshift = atof(argv[1]);
	nx = atoi(argv[2]);

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
	MPI_Status status;

	MPI_Barrier(MPI_COMM_WORLD);
	printf("P%d passed\n",myid);


	sprintf(infile,"cplMergingTree.dat");
	MajorGid *smbp;
	IDTYPE npsmbp = readmergerdata(infile, &smbp );

	IDTYPE gmin, gmax, lgmin, lgmax, flgmin, flgmax;
	GetMinMaxofGid(smbp, npsmbp, &gmin, &gmax);

	GetLocalMinMaxofGid(smbp, npsmbp, &flgmin, &flgmax, gmin, gmax);


	IDTYPE icount = 0;
	for(i=0;i<npsmbp;i++) {
		if(smbp[i].majorgid[iredshift] <0) {
			smbp[i].mgid = LONG_MAX;
			icount ++;
		}
		else smbp[i].mgid = smbp[i].majorgid[iredshift];
	}
	printf("P%d has negative halo value %ld\n",myid, icount);


	AdGetLocalMinMaxofGid(smbp, npsmbp, &lgmin, &lgmax, gmin, gmax);

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
		AdGetLocalMinMaxofGid(smbp, npsmbp, &lgmin, &lgmax, gmin, gmax);
		MyFixedRangePqsort(smbp,npsmbp,MajorGid,IDTYPE,mgid,lgmin,lgmax,MPI_COMM_WORLD);
		MyFixedRangePqsort(newgid, npnewgid, GID, IDTYPE, gid, lgmin, lgmax, MPI_COMM_WORLD);
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


	MyFixedRangePqsort(smbp,npsmbp,MajorGid,IDTYPE,gid,flgmin,flgmax,MPI_COMM_WORLD);
	
	sprintf(outfile,"SortByMajorGidMergingTree.%.5d.dat", iredshift);
	writehaloQ(smbp, npsmbp, infile, outfile, flgmin,flgmax);


	if(myid==0) printf("Exit with success\n");
	MPI_Finalize();
	return 0;

}

