/* This code outputs the major merging tree of the most massive halos at the input time step */
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
#include "pqsort.h"
#include "addQ.h"

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

double Msun=1.989E33L;
double G= 6.672E-8L;
double pc= 3.08567802E18L;
double PI= 3.14159265354979L;
double onesolarmass=1.989E33L;

int nsteps[NSTEPS] ={
	167 ,191 ,213 ,235 ,256 ,276 ,295 ,314 ,332 ,350 ,368 ,385 ,401 ,418 ,434 ,450 ,466 ,481 ,496
		,511 ,526 ,541 ,556 ,570 ,584 ,598 ,612 ,626 ,640 ,654 ,667 ,681 ,694 ,708 ,721 ,734 ,747 ,760 ,773
		,786 ,799 ,812 ,825 ,837 ,850 ,863 ,875 ,888 ,900 ,913 ,925 ,938 ,950 ,963 ,975 ,987 ,1000 ,1012 ,1024
		,1037 ,1049 ,1061 ,1074 ,1086 ,1098 ,1110 ,1123 ,1135 ,1147 ,1160 ,1172 ,1184 ,1197 ,1209 ,1221 ,1234 ,1246 ,1258 ,1271
		,1283 ,1296 ,1308 ,1321 ,1333 ,1358 ,1371 ,1383 ,1396 ,1409 ,1421 ,1434 ,1447 ,1460 ,1472 ,1485 ,1498 ,1511 ,1524 ,1537
		,1550 ,1563 ,1576 ,1590 ,1603 ,1616 ,1629 ,1656 ,1669 ,1683 ,1696 ,1710 ,1737 ,1751 ,1765 ,1779 ,1793 ,1807 ,1821 ,1835
		,1849 ,1863 ,1877 ,1892 ,1906 ,1920 ,1935 ,1949 ,1964 ,1979 ,1981
};

/*
int nsteps[NSTEPS] ={
	167 ,191 ,213 ,235 ,256 ,276 ,295 ,314 ,332 ,350 ,368 ,385 ,401 ,418 ,434 ,450 
		,466 ,481 ,496 ,511 ,526 ,541 ,556 ,570 ,584 ,598 ,612 ,626 ,640 ,654 ,667 ,681 ,694 
		,708 ,721 ,734 ,747 ,760 ,773 ,786 ,799 ,812 ,825 ,837 ,850 ,863 ,875 ,888 ,900 ,913
		,925 ,938 ,950 ,963 ,975 ,987 ,1000 ,1012 ,1024 ,1037 ,1049 ,1061 ,1074 ,1086 ,1098 ,1110 ,1123 
		,1135 ,1147 ,1160 ,1172 ,1184 ,1197 ,1209 ,1221 ,1234 ,1246 ,1258 ,1271 ,1283 ,1296 ,1308 ,1321 ,1333
		,1346 ,1358 ,1371 ,1383 ,1396 ,1409 ,1421 ,1434 ,1447 ,1460 ,1472 ,1485 ,1498 ,1511 ,1524 ,1537 ,1550
		,1563 ,1576 ,1590 ,1603 ,1616 ,1629 ,1643 ,1656 ,1669 ,1683 ,1696 ,1710 ,1724 ,1737 ,1751 ,1765 ,1779
		,1793 ,1807 ,1821 ,1835 ,1849 ,1863 ,1877 ,1892 ,1906 ,1920 ,1935 ,1949 ,1964 ,1979 ,1981
};
*/




int myid, nid;
double pmass,minmass,rnx;
float boxsize,hubble,npower,omep,omepb,omeplam,bias,amax,astep,anow, wlam0, wlam1;
int nx,nspace;


size_t initshift;


#define MIN(a,b) ((a)<(b) ? (a): (b))
#define MAX(a,b) ((a)>(b) ? (a): (b))
#define pow2(A) ((A)*(A))
#define pow3(A) ((A)*(A)*(A))
#define sqr(A) ((A)*(A))



typedef struct MergerAddQ{
	float haloq[NSTEPS];
} MergerAddQ;

IDTYPE readhalobindata(char *infile, HID **hid){
	float iredshift;
	IDTYPE tnp,np,i,j,k;
	IDTYPE stnp, wtnp, ftnp;
	size_t nshift = 0;
	if(myid==0){
		FILE *fp = fopen(infile,"r");
		fseek(fp, 0, SEEK_END);
		long len = (long) ftell(fp);
		fclose(fp);
		len = len - nshift;
		tnp = len/sizeof(HaloProp);
		printf("Total number of Halo-Quantity data is %ld ::: %ld\n", tnp, len%(long)sizeof(HaloProp));
	}
	MPI_Bcast(&tnp, sizeof(IDTYPE), MPI_BYTE, 0, MPI_COMM_WORLD);
	wtnp = (tnp + nid-1)/nid;
	stnp = myid*wtnp;
	ftnp = (myid+1)*wtnp;
	ftnp = MIN(ftnp, tnp);

	wtnp = ftnp - stnp;
	IDTYPE isend, iget ;
	isend = iget =0;
	int WGroupSize = 10000000;/* Mandatory */
	int src = myid-1;
	int tgt = myid+1;
	int itag= 0;
	MPI_Status status;


	HID *nbp = (HID*)malloc(sizeof(HID) * wtnp);
	*hid = nbp;
	printf("P%d has min mass %g\n", myid, minmass);

	if(myid==-1){
		int kk=1;
		while(kk){
			i = i;
		}
	}

	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_LONG,src,itag,MPI_COMM_WORLD,&status);
	FILE *fp = fopen(infile,"r");
	nshift += sizeof(HaloProp)*stnp;

	fseek(fp, nshift, SEEK_SET);
	HaloProp *bp = (HaloProp*)malloc(sizeof(HaloProp)*100000);
	IDTYPE itnp = 0;
	IDTYPE icount = 0;
	if(wtnp >0){
		while((np = fread(bp,sizeof(HaloProp),100000,fp))){
			for(i=0;i<np;i++){
				icount ++;
				if(bp[i].mass >= minmass) {
					nbp[itnp].hid = iget;
					nbp[itnp].sigv = bp[i].sigv;
					itnp ++;
					iget ++;
				}
				if(icount >= wtnp) goto out;
			}
		}
	}
	else{
		*hid = NULL;
	}
out:
	fclose(fp);
	free(bp);

	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&iget,1,MPI_LONG,tgt,itag,MPI_COMM_WORLD);

	return itnp;
}
typedef struct MergeHistory{
    SavedHaloType step[NSTEPS];
} MergeHistory;



IDTYPE readmergerdata(char *infile, MBP **mbp){
	float iredshift;
	IDTYPE tnp,np,i,j,k;
	IDTYPE stnp, wtnp;
	size_t nshift = sizeof(float)*10 + sizeof(int);
	nshift = sizeof(float)*(10+NSTEPS) + sizeof(int)*(1+NSTEPS);

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
		fread(&wlam0,sizeof(float),1,fp);
		fread(&wlam1,sizeof(float),1,fp);
		fread(&bias,sizeof(float),1,fp);
		nspace = 1;

		rnx = nx;


		pmass = 3./8./M_PI/G*sqr(100.E5*hubble)*1.E6*(pc/onesolarmass);
		pmass = pmass*pow3(boxsize/hubble)/pow3(rnx/nspace)*omep;

		minmass = 30.L*pmass * 0.999L;
		printf("The minimum halo mass is %g :: %g %d %g\n",minmass, boxsize, nx, omep);

		fclose(fp);
		len = len - nshift;
		tnp = len/(sizeof(MergeHistory));
		printf("Total number of MBP is %ld :: %ld\n", tnp, len%(sizeof(MergeHistory)));
	}
	MPI_Bcast(&tnp, sizeof(IDTYPE), MPI_BYTE, 0, MPI_COMM_WORLD);
	MPI_Bcast(&minmass, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	wtnp = (tnp + nid-1)/nid;
	stnp = myid*wtnp;

	wtnp = MIN(wtnp, tnp - stnp);
	int isend, iget ;
	isend = iget =1;
	int WGroupSize = WGROUPSIZE;
	int src = myid-1;
	int tgt = myid+1;
	int itag= 0;
	MPI_Status status;


	MBP *nbp = (MBP*)malloc(sizeof(MBP) * wtnp);
	*mbp =  nbp;


	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	FILE *fp = fopen(infile,"r");
	nshift += sizeof(MergeHistory)*stnp;
	printf("P%d tries to read merger data :::: %ld\n",myid, wtnp);

	fseek(fp, nshift, SEEK_SET);
	MergeHistory *bp = (MergeHistory*)malloc(sizeof(MergeHistory)*1000);
	IDTYPE itnp = 0;
	IDTYPE icount = 0;
	while((np = fread(bp,sizeof(MergeHistory),1000,fp))){
		for(i=0;i<np;i++){
			nbp[itnp].gid = stnp + itnp;
			for(j=0;j<NSTEPS;j++){
				if(bp[i].step[j].mbp<=0){
					nbp[itnp].hid[j] = -1;
				}
				else{
					nbp[itnp].hid[j] = bp[i].step[j].nowhid;
					if(nbp[itnp].hid[j] ==0) {
						icount ++;
					}
				}
			}
			itnp ++;
			if(itnp >= wtnp) goto out;
		}
	}
out:
	printf("P%d Total number of Icount = %ld\n",myid, icount);
	fclose(fp);
	free(bp);

	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);


	return wtnp;
}

void writehaloQ(MergerAddQ *haloq, char *outfile, IDTYPE gmin, IDTYPE gmax){ 
	IDTYPE i,j,k; 
	int isend, iget ; 
	isend = iget =1; 
	int WGroupSize = WGROUPSIZE; 
	int src = myid-1; 
	int tgt = myid+1; 
	int itag = 18; 
	MPI_Status status; 
	if(myid==nid-1){
		FILE *fp = fopen(outfile, "w");
		long shift = sizeof(MergerAddQ)*gmin;
		printf("P%d is writing %s with shift %ld\n",myid, outfile, shift);
		fseek(fp, shift, SEEK_SET);
		fwrite(haloq, sizeof(MergerAddQ),(gmax-gmin), fp);
		fclose(fp);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status); 
	if(myid!=nid-1){
		FILE *fp = fopen(outfile, "r+");
		long shift = sizeof(MergerAddQ)*gmin;
		fseek(fp, shift, SEEK_SET);
		fwrite(haloq, sizeof(MergerAddQ),(gmax-gmin), fp);
		fclose(fp);
	}
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
}

void GetMinMaxofGid(MBP *hid,IDTYPE npnow, IDTYPE *gmin, IDTYPE *gmax){
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
void GetMinMaxofHid(HID *hid,IDTYPE npnow, IDTYPE *hmin, IDTYPE *hmax){
	IDTYPE i,j,k;
	IDTYPE rmin,rmax,tmin,tmax;
	rmin = INT_MAX;
	rmax = INT_MIN;
	for(i=0;i<npnow;i++){
		if(rmin > hid[i].hid) rmin = hid[i].hid;
		if(rmax < hid[i].hid) rmax = hid[i].hid;
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
	*hmin = rmin;
	*hmax = rmax;
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


	nx = atof(argv[1]);

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
	float nredshift,iredshift;
	IDTYPE npmbp;
	MBP *mbp;
	MPI_Status status;


	if(myid==-1){
		int kk=1;
		while(kk){
			i = i;
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);


	sprintf(infile,"MergingTree.dat");

	npmbp = readmergerdata(infile, &mbp );
	IDTYPE min,max,gmin,gmax;
	GetMinMaxofGid(mbp, npmbp, &gmin,&gmax);
	min = gmin + myid*( (gmax-gmin-1)/nid +1);
	max = gmin + (myid+1)*((gmax-gmin-1)/nid+1);
	min = MIN(min,gmax);
	max = MIN(max,gmax);

	MergerAddQ *haloq = (MergerAddQ*)calloc((max-min+1), sizeof(MergerAddQ));


	if(myid == nid-1) max = gmax + 1;


	for(i=0; i<NSTEPS  ;i++){

		IDTYPE inp;
		sprintf(infile,"cplFoFHaloQuantities%.5d.bin",nsteps[i]);
		printf("Opening %s\n",infile);
		HID *hid;
		IDTYPE nphid, hmin,hmax,lhmin,lhmax;
		nphid= readhalobindata(infile,&hid );
		GetMinMaxofHid(hid, nphid, &hmin, &hmax);
		lhmin = hmin + myid*( (hmax-hmin -1)/nid + 1);
		lhmax = hmin + (myid+1)*( (hmax-hmin -1)/nid + 1);
		lhmin = MIN(lhmin, hmax);
		lhmax = MIN(lhmax, hmax);
		if(myid == nid-1) lhmax = hmax + 1;


		IDTYPE npwmbp = 0;
		for(j=0;j<npmbp;j++) if(mbp[j].hid[i] >=0) npwmbp ++;
		HID *wmbp = (HID*)malloc(sizeof(HID)*npwmbp);
		npwmbp = 0;
		for(j=0;j<npmbp;j++) if(mbp[j].hid[i] >=0) {
			wmbp[npwmbp].hid = mbp[j].hid[i];
			wmbp[npwmbp].gid = mbp[j].gid;
			npwmbp ++;
		}
		IDTYPE whmin, whmax;
		GetMinMaxofHid(wmbp, npwmbp, &whmin, &whmax);
		IDTYPE tnphid, tnpwmbp;
		MPI_Reduce(&nphid,&tnphid, 1, MPI_LONG,MPI_SUM,0, MPI_COMM_WORLD);
		MPI_Reduce(&npwmbp,&tnpwmbp, 1, MPI_LONG,MPI_SUM,0, MPI_COMM_WORLD);


		if(myid==0) {
			printf("##############################################\n");
			printf("TO CHECK DIFFERENCE %ld - %ld : %ld :::: %ld - %ld : %ld\n", hmin, hmax, tnphid,
					whmin, whmax, tnpwmbp);
			printf("##############################################\n");
		}
		debug(wmbp, npwmbp, 0);
		MyFixedRangePqsort(wmbp,npwmbp,HID,IDTYPE,hid,lhmin,lhmax,MPI_COMM_WORLD);
		debug(wmbp, npwmbp, 1);
		MyFixedRangePqsort(hid,nphid,HID,IDTYPE,hid,lhmin,lhmax, MPI_COMM_WORLD);
		debug(wmbp, npwmbp, 2);
		j = 0;
		for(k=0;k<npwmbp;k++){
			for(;j<nphid;j++){
				if(wmbp[k].hid == hid[j].hid) {
					wmbp[k].sigv = hid[j].sigv;
					break;
				}
			}
		}
		debug(wmbp, npwmbp, 3);

		MyFixedRangePqsort(wmbp,npwmbp,HID,IDTYPE,gid,min,max,  MPI_COMM_WORLD);
		debug(wmbp, npwmbp, 4);
		for(j=0;j<npwmbp;j++){
			size_t ishift = wmbp[j].gid - min;
			haloq[ishift].haloq[i] = wmbp[j].sigv;
		}

		free(hid);
		free(wmbp);
		printf("P%d End of %ld :: %d\n",myid, i, NSTEPS);
	}

	sprintf(outfile,"MergerSgiv.dat");
	writehaloQ(haloq, outfile, min,max);


	MPI_Finalize();
	return 0;

}

