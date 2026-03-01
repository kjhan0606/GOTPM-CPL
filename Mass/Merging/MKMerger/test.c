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

#define GROUPID(a,b) ((a)/(b))
#define RANKINGROUP(a,b) ((a)%(b))

double Msun=1.989E33L;
double G= 6.672E-8L;
double pc= 3.08567802E18L;
double PI= 3.14159265354979L;



#define NSTEPS 131
int nsteps[NSTEPS] ={
167 , 191 , 213 , 235 , 256 , 276 , 295 , 314 , 332 , 350 , 368 , 418 , 434 , 450 , 466 , 481 , 496 , 511 , 526 , 541 , 556 , 570 , 584 , 598 , 612 , 626 ,
640 , 654 , 667 , 681 , 694 , 708 , 721 , 734 , 747 , 760 , 773 , 786 , 799 , 812 , 825 , 837 , 850 , 863 , 875 , 888 , 900 , 913 , 925 , 938 , 950 , 963 , 975 , 987 , 1000,
1012, 1024, 1037, 1049, 1061, 1074, 1086, 1098, 1110, 1123, 1135, 1147, 1160, 1172, 1184, 1197, 1209, 1221, 1234, 1246, 1258, 1271, 1283, 1296,
1308, 1321, 1333, 1346, 1358, 1371, 1383, 1396, 1409, 1421, 1434, 1447, 1460, 1472, 1485, 1498, 1511, 1524, 1537, 1550, 1563, 1576, 1590, 1603, 1616, 1629, 1643, 1656,
1669, 1683, 1696, 1710, 1724, 1737, 1751, 1765, 1779, 1793, 1807, 1821, 1835, 1849, 1863, 1877, 1892, 1906, 1920, 1935, 1949, 1964, 1979, 1981
};
/*
#define NSTEPS 140
int nsteps[NSTEPS] ={
176,198,218,237,256,274,291,308,325,341,357,373,388,403,418,433,447,461,475,489,503,
517,530,543,556,570,583,595,608,621,634,646,659,671,683,695,708,720,732,744,756,768,780,792,
804,815,827,839,851,862,874,886,897,909,921,932,944,955,967,979,990,1002,1013,1025,1036,1048,1060,
1071,1083,1094,1106,1118,1129,1141,1153,1164,1176,1188,1200,1211,1223,1235,1247,1259,1271,1282,1294,1306,1318,1330,
1343,1355,1367,1379,1391,1403,1416,1428,1440,1453,1465,1478,1490,1503,1516,1528,1541,1554,1567,1579,1592,1605,1618,
1632,1645,1658,1671,1684,1698,1711,1725,1738,1752,1765,1779,1793,1807,1821,1835,1849,1863,1877,1891,1906,1920,1934,
1949,1964,1978,1981};
*/

/*
#define NSTEPS 126
int nsteps[NSTEPS] ={
178,  203,  227,  250,  272,  293,  314,  334,  353,  372,  391,
409,  427,  444,  461,  478,  495,  511,  527,  543,  559,  574,  590,
605,  620,  635,  650,  665,  679,  694,  708,  722,  737,  751,  765,
779,  792,  806,  820,  833,  847,  861,  874,  887,  901,  914,  927,
941,  954,  967,  980,  993,  1006,  1019,  1032,  1045,  1058,  1071,  1084,
1097,  1110,  1123,  1136,  1149,  1162,  1175,  1187,  1200,  1213,  1226,  1239,
1252,  1265,  1278,  1291,  1304,  1316,  1329,  1342,  1355,  1368,  1381,  1394,
1407,  1420,  1434,  1447,  1460,  1473,  1486,  1499,  1513,  1526,  1539,  1552,
1566,  1579,  1593,  1606,  1620,  1633,  1647,  1660, 1674, 1688, 1701, 1715,
1729, 1743, 1757, 1771, 1785, 1799, 1813, 1827, 1841, 1855, 1869, 1884, 
1898, 1913, 1927, 1942, 1956, 1971, 1981
};
*/


int myid, nid;
double pmass;
float boxsize,hubble,npower,omep,omepb,omeplam,bias,amax,astep,anow;
int nx,nspace;


size_t initshift;


#define MIN(a,b) ((a)<(b) ? (a): (b))
#define MAX(a,b) ((a)>(b) ? (a): (b))



IDTYPE readdata(char *infile, char *infile2, char *infile3, char *infile4,SimplifiedHalo **simhalo ){
	float iredshift;
	IDTYPE tnp,np,i,j,k;
	FILE *fp = fopen(infile,"r");
	if(myid==0){
		fread(&iredshift, sizeof(float),1,fp);
		fread(&tnp,sizeof(IDTYPE),1,fp);
		printf("tnp = %ld in %s\n",tnp, infile);
	}
	MPI_Bcast(&tnp,sizeof(IDTYPE), MPI_BYTE,0, MPI_COMM_WORLD);
	IDTYPE nstep = (tnp-1)/nid + 1;
	IDTYPE nstart = MIN(tnp,nstep*myid);
	np =  MIN(tnp, nstep*(myid+1)) - MIN(tnp, nstep*myid);

	size_t nshift = sizeof(float)+sizeof(IDTYPE) + sizeof(TrHalo)*nstart;


	int isend, iget ;
	isend = iget =1;
	int WGroupSize = WGROUPSIZE;
	int src = myid-1;
	int tgt = myid+1;
	int itag= 0;
	MPI_Status status;


	void readsiminfo( char *);
	readsiminfo(infile3);



	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);
	fseek(fp, nshift, SEEK_SET);
	TrHalo *bp = (TrHalo*)malloc(sizeof(TrHalo)*np);
	SimplifiedHalo *shalo = NULL;
	if(np) shalo = (SimplifiedHalo*)malloc(sizeof(SimplifiedHalo)*np);
	
	/*
	printf("P%d has np = %ld : shalo = %p\n",myid,np,shalo);
	*/
	fread(bp,sizeof(TrHalo),np,fp);
#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		/*
		shalo[i].nowaid = nstart + i;
		*/
		shalo[i].nowaid = bp[i].aid;
		shalo[i].downaid = bp[i].downaid;
		shalo[i].upaid = bp[i].upaid;
		shalo[i].basichalo.mbp = bp[i].mbp;
		shalo[i].globalid = -1;
		shalo[i].basichalo.mbp = bp[i].mbp;
		shalo[i].basichalo.nowhid = (shalo[i].nowhid = bp[i].nowhid);
		shalo[i].basichalo.mass = bp[i].mass*hubble;
		shalo[i].basichalo.fofhmass = bp[i].fofhmass*hubble;
		shalo[i].basichalo.statusflag = bp[i].statusflag;
		shalo[i].basichalo.linkflag = bp[i].linkflag;

	}
	free(bp);
	fclose(fp);
	MbpRV *rv = (MbpRV*)malloc(sizeof(MbpRV)*np);
	MbpRV *hrv = (MbpRV*)malloc(sizeof(MbpRV)*np);
	fp = fopen(infile2,"r");
	nshift = sizeof(MbpRV)*nstart;
	fseek(fp,nshift,SEEK_SET);
	fread(rv,sizeof(MbpRV),np,fp);
	fclose(fp);

	fp = fopen(infile4,"r");
	nshift = sizeof(MbpRV)*nstart;
	fseek(fp,nshift,SEEK_SET);
	fread(hrv,sizeof(MbpRV),np,fp);
	fclose(fp);



	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid)
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);



#ifdef _OPENMP
#pragma omp parallel for
#endif
	for(i=0;i<np;i++){
		shalo[i].basichalo.x = rv[i].x;
		shalo[i].basichalo.y = rv[i].y;
		shalo[i].basichalo.z = rv[i].z;
		shalo[i].basichalo.vx = rv[i].vx;
		shalo[i].basichalo.vy = rv[i].vy;
		shalo[i].basichalo.vz = rv[i].vz;
		shalo[i].basichalo.hx = hrv[i].x;
		shalo[i].basichalo.hy = hrv[i].y;
		shalo[i].basichalo.hz = hrv[i].z;
		shalo[i].basichalo.hvx = hrv[i].vx;
		shalo[i].basichalo.hvy = hrv[i].vy;
		shalo[i].basichalo.hvz = hrv[i].vz;
		if(shalo[i].basichalo.mbp == 7250204869L)
			printf("-P%d has found  %ld : %g : %s\n",myid,shalo[i].basichalo.mbp, shalo[i].basichalo.x, infile);
	}
	free(rv);





	/*
	HidMass *hidmass;
	IDTYPE readhidmass( char *, HidMass **);
	IDTYPE nnowhid = readhidmass(infile3,&hidmass);

	IDTYPE min,max,tmin,tmax;
	min = INT_MAX;
	max = INT_MIN;
	for(i=0;i<nnowhid;i++){
		if(hidmass[i].hid < min) min = hidmass[i].hid;
		if(hidmass[i].hid > max) max = hidmass[i].hid;
	}
	MPI_Reduce(&min,&tmin,1,MPI_IDTYPE,MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&max,&tmax,1,MPI_IDTYPE,MPI_MAX, 0, MPI_COMM_WORLD);
	for(i=0;i<np;i++){
		if(shalo[i].nowhid < min) min = shalo[i].nowhid;
		if(shalo[i].nowhid > max) max = shalo[i].nowhid;
	}
	MPI_Reduce(&min,&tmin,1,MPI_IDTYPE,MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&max,&tmax,1,MPI_IDTYPE,MPI_MAX, 0, MPI_COMM_WORLD);
	if(myid==0){
		min = tmin;
		max = tmax;
	} 
	MPI_Bcast(&min,1,MPI_IDTYPE,0,MPI_COMM_WORLD);
	MPI_Bcast(&max,1,MPI_IDTYPE,0,MPI_COMM_WORLD);
	IDTYPE gmin,gmax;
	gmin = min;
	gmax = max;
	IDTYPE step = (gmax-gmin-1)/nid + 1;
	min = MIN(gmax, gmin + myid*step);
	max = MIN(gmax, gmin + (myid+1)*step);
	if(myid == nid-1) max = gmax + 1;




	MyFixedRangePqsort(shalo,np,SimplifiedHalo,IDTYPE,nowhid,min,max,
			MPI_COMM_WORLD);
	MyFixedRangePqsort(hidmass,nnowhid,HidMass,IDTYPE,hid,min,max,
			MPI_COMM_WORLD);
	i = 0;
	for(j=0;j<nnowhid;j++){
		for(;i<np;i++){
			if(shalo[i].nowhid == hidmass[j].hid){
				shalo[i].basichalo.fofhmass = hidmass[j].mass;
				if(shalo[i].basichalo.mass > 1.5E13 &&  hidmass[j].mass< 2E12){
					char aaa;
					if( IS_MAJOR_SATELLITE(shalo,i,basichalo.statusflag)) aaa = 'M';
					else if( IS_MINOR_SATELLITE(shalo,i,basichalo.statusflag)) aaa = 'S';
					printf("P%d has strange value %c %g %g\n",myid, aaa,shalo[i].basichalo.mass, hidmass[j].mass);
					exit(19);
				}
			}
			else if(shalo[i].nowhid > hidmass[j].hid) break;
		}
	}
	free(hidmass);
	*/







	*simhalo = shalo;


	return np;
}
void readsiminfo(char *infile2){
        IDTYPE i,j,k;
        IDTYPE maxhidmass = 10000000,nhidmass=0,np;
        int isend, iget ;
        isend = iget =1;
        int WGroupSize = WGROUPSIZE;
        int src = myid-1;
        int tgt = myid+1;
        int itag = 10;
        MPI_Status status;
        if(RANKINGROUP(myid,WGroupSize) != 0 )
                MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);


        FILE *fp = fopen(infile2, "r");
        {
                FILE *fp = fopen(infile2,"r");
                fread(&boxsize,sizeof(float),1,fp);
                fread(&hubble,sizeof(float),1,fp);
                fread(&npower,sizeof(float),1,fp);
                fread(&omep,sizeof(float),1,fp);
                fread(&omepb,sizeof(float),1,fp);
                fread(&omeplam,sizeof(float),1,fp);
                fread(&bias,sizeof(float),1,fp);
                fread(&nx,sizeof(int),1,fp);
                fread(&nspace,sizeof(int),1,fp);
                fread(&amax,sizeof(float),1,fp);
                fread(&astep,sizeof(float),1,fp);
                fread(&anow,sizeof(float),1,fp);
                pmass = 3.L/8.L/PI/G*pow(100.E5,2.L)*1.E6*pc;
                pmass = pmass*pow((double)boxsize,3.L)/pow(nx/nspace,3.L)*omep/Msun;
	}
	fclose(fp);
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
}



IDTYPE readhidmass(char *infile2,HidMass **hidmass){
	IDTYPE i,j,k;
	IDTYPE maxhidmass = 10000000,nhidmass=0,np;
	if(myid==0){
		struct stat st;
		stat(infile2,&st);
		np = (st.st_size-sizeof(float)*10-sizeof(int)*2)/sizeof(Halo);
		printf("%s : np = %ld : %ld %ld\n",infile2,np, st.st_size, sizeof(Halo));
	}
	MPI_Bcast(&np,sizeof(IDTYPE),MPI_BYTE,0,MPI_COMM_WORLD);



	IDTYPE nstep = (np-1)/nid +1;

	IDTYPE mp = MIN(np,nstep*(myid+1)) - MIN(np,nstep*myid);
	IDTYPE nstart = MIN(np,nstep*myid);


	int isend, iget ; 
	isend = iget =1; 
	int WGroupSize = WGROUPSIZE; 
	int src = myid-1; 
	int tgt = myid+1; 
	int itag = 10;
	MPI_Status status;
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);


	FILE *fp = fopen(infile2, "r");
	{
		FILE *fp = fopen(infile2,"r");
		int nshift = sizeof(float)*10 + sizeof(int)*2 + nstart*sizeof(Halo);
		fread(&boxsize,sizeof(float),1,fp);
		fread(&hubble,sizeof(float),1,fp);
		fread(&npower,sizeof(float),1,fp);
		fread(&omep,sizeof(float),1,fp);
		fread(&omepb,sizeof(float),1,fp);
		fread(&omeplam,sizeof(float),1,fp);
		fread(&bias,sizeof(float),1,fp);
		fread(&nx,sizeof(int),1,fp);
		fread(&nspace,sizeof(int),1,fp);
		fread(&amax,sizeof(float),1,fp);
		fread(&astep,sizeof(float),1,fp);
		fread(&anow,sizeof(float),1,fp);
		pmass = 3.L/8.L/PI/G*pow(100.E5,2.L)*1.E6*pc;
		pmass = pmass*pow((double)boxsize,3.L)/pow(nx/nspace,3.L)*omep/Msun;

		fseek(fp,nshift,SEEK_SET);
		Halo *halo = (Halo*)malloc(sizeof(Halo)*mp);
		IDTYPE mmp = fread(halo,sizeof(Halo),mp,fp);
		fclose(fp);
		if(mmp != mp){
			printf("P%d has reading error %ld : %ld\n",myid,mp,mmp);
			exit(99);
		}
		*hidmass = (HidMass *)malloc(sizeof(HidMass)*mp);
		j = 0;
		for(i=0;i<mp;i++){
			if(halo[i].np >= NMINMEMBERS) {
				/*
				(*hidmass)[j].hid = nstart*myid + j;
				*/
				(*hidmass)[j].mass = halo[i].np*pmass;
				j++;
			}
			/*
			printf("P%d has %ld halo.np -> %ld\n",myid, i, halo[i].np);
			*/
		}
		printf("P%d read fof mass data: %ld %ld->%ld : tnp %ld\n",myid,nstart, mp,j, np);
		mp = j;
		free(halo);
	}

	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);

	IDTYPE globalhid = 0;
	if(myid != 0) MPI_Recv(&globalhid, 1, MPI_IDTYPE, (myid-1), 1, MPI_COMM_WORLD, &status);

	for(i=0;i<mp;i++) {
		(*hidmass)[i].hid = globalhid++;
	}
	if(myid != nid-1) MPI_Send(&globalhid, 1, MPI_IDTYPE, myid+1, 1, MPI_COMM_WORLD);



	return mp;
}
void determinemajorglobalid(SimplifiedHalo *shnow, IDTYPE npnow){
	IDTYPE i,j,k;

	for(i=0;i<npnow;){
		for(j=i+1;j<npnow;j++){
			if(shnow[j].nowhid != shnow[i].nowhid) break;
		}
		float maxmass = shnow[i].basichalo.mass;
		IDTYPE imax=i;
		for(k=i+1;k<j;k++){
			if(shnow[k].basichalo.mass > maxmass){
				imax = k;
				maxmass = shnow[k].basichalo.mass;
			}
		}
		for(k=i;k<j;k++){
			shnow[k].basichalo.majorglobalid = shnow[imax].globalid;
		}
		i = j;
	}
}

IDTYPE globalid = 0;
void getminmaxofupnowaid(SimplifiedHalo *shnow,IDTYPE npnow, SimplifiedHalo *shup,IDTYPE npup, IDTYPE *min, IDTYPE *max){
	IDTYPE i,j,k;
	IDTYPE rmin,rmax,tmin,tmax;
	rmin = INT_MAX;
	rmax = INT_MIN;
	for(i=0;i<npnow;i++){
		if(rmin > shnow[i].upaid) rmin = shnow[i].upaid;
		if(rmax < shnow[i].upaid) rmax = shnow[i].upaid;
	}
	for(i=0;i<npup;i++){
		if(rmin > shup[i].nowaid) rmin = shup[i].nowaid;
		if(rmax < shup[i].nowaid) rmax = shup[i].nowaid;
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
	*min = rmin;
	*max = rmax;
}


/*
void MergeFiles(char *outfile,IDTYPE globalid){
	IDTYPE i,j,k,gidblock=400000;

	FILE *wp = fopen(outfile,"a");

	SavedHaloType *bp = (SavedHaloType*)malloc(sizeof(SavedHaloType)*NSTEPS*gidblock);
	SavedHaloType *btmp = (SavedHaloType*)malloc(sizeof(SavedHaloType)*gidblock);
	for(j=0;j<globalid;j+= gidblock){
		IDTYPE readblock = MIN(gidblock, globalid-j);
		for(i=0;i<NSTEPS;i++){
			char out1[190];
			sprintf(out1,"%s%.5d",outfile,(int)i);
			FILE *fp = fopen(out1,"r");
			printf("P%d is reading %s with gid= %ld\n",myid, out1,j);
			fseek(fp,j*sizeof(SavedHaloType),SEEK_SET);
			long long mmp = fread(btmp,sizeof(SavedHaloType),readblock,fp);
			for(k=0;k<mmp;k++){
				bp[i + NSTEPS*k] = btmp[k];
			}
			fclose(fp);
		}
		fwrite(bp,sizeof(SavedHaloType),NSTEPS*gidblock,wp);
	}
	fclose(wp);
	free(bp);
	free(btmp);
	
}
*/


void dumpheader(char *outfile){
	FILE *wp;

	/*
	wp = fopen(outfile,"w");
	fwrite(&boxsize,sizeof(float),1,wp);
	fwrite(&hubble,sizeof(float),1,wp);
	fwrite(&npower,sizeof(float),1,wp);
	fwrite(&omep,sizeof(float),1,wp);
	fwrite(&omepb,sizeof(float),1,wp);
	fwrite(&omeplam,sizeof(float),1,wp);
	fwrite(&bias,sizeof(float),1,wp);
	fwrite(&astep,sizeof(float),1,wp);
	int msteps = NSTEPS;
	fwrite(&msteps,sizeof(int),1,wp);
	fwrite(nsteps,sizeof(int),msteps,wp);
	int i;
	for(i=0;i<msteps;i++){
		float red = amax/(1+astep*(nsteps[i]-1)) - 1;
		fwrite(&red,sizeof(float),1,wp);
	}
	fclose(wp);
	*/
}
void debug(SimplifiedHalo *arr, IDTYPE np, char *infile){
	IDTYPE i=0;
	for(i=0;i<np;i++){
		if(arr[i].basichalo.mbp == 7250204869L){
			printf("+P%d has found  %ld : %g : %s\n",myid,arr[i].basichalo.mbp, arr[i].basichalo.x, infile);
			break;
		}
	}
	fflush(stdout);
}

int main(int argc, char **argv){
	IDTYPE i,j,k;
	IDTYPE np1,np2,ip;
	TrHalo *nmbp,imbp;
	IDTYPE bmem, amem;
	MbpRV imbprv;
	FILE *fp1;
	HidMass *hidmass;
	char outfile[190];


	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	if(myid==0){
		printf("###########################################################################\n");
		printf("This program generates galaxies from the merging tree.\n");
		printf("###########################################################################\n");
	}
        {
                char hostname[190];
                hostname[189] = '\0';
                gethostname(hostname,189);
                printf("P%d host name is %s with pid= %d\n",myid,hostname, getpid());
        }



	sprintf(outfile,"MergingTree.dat");

	char infile[190],infile1[190],infile2[190],infile3[100];
	FILE *fp;
	float nredshift,iredshift;
	IDTYPE npup,npnow;
	SimplifiedHalo *shup,*shnow;
	MPI_Status status;



	/*
	if(myid==9){
		int kk=1;
		while(kk){
			i = i;
		}
	}
	*/


	sprintf(infile,"HaloLinkedList.%.5d",nsteps[0]);
	sprintf(infile1,"HaloRVList.%.5d",nsteps[0]);
	sprintf(infile2,"FoF_halo_cat.%.5d",nsteps[0]);
	sprintf(infile3,"HaloRVList_halo.%.5d",nsteps[0]);
	npnow = readdata(infile, infile1, infile2,infile3,&shnow );
	debug(shnow,npnow, infile);

	/* Setting initial global id to the first slice data */
	if(myid!=0) MPI_Recv(&globalid,sizeof(IDTYPE),MPI_BYTE,myid-1,0,MPI_COMM_WORLD,&status);
	for(i=0;i<npnow;i++){
		shnow[i].globalid = globalid;
		globalid++;
	}
	if(myid != nid-1) MPI_Send(&globalid,sizeof(IDTYPE),MPI_BYTE,myid+1,0,MPI_COMM_WORLD);
	MPI_Bcast(&globalid,sizeof(IDTYPE),MPI_BYTE,nid-1,MPI_COMM_WORLD);


	determinemajorglobalid(shnow, npnow);




	initshift = sizeof(float)*8 + sizeof(int)*1 + (sizeof(int)+sizeof(float))*NSTEPS;





	/*
	if(myid==1){
		int kk=1;
		while(kk){
			i = i;
		}
	}
	*/


	for(i=1; i<NSTEPS  ;i++){
		IDTYPE inp;
		sprintf(infile,"HaloLinkedList.%.5d",nsteps[i]);
		sprintf(infile1,"HaloRVList.%.5d",nsteps[i]);
		sprintf(infile2,"FoF_halo_cat.%.5d",nsteps[i]); 
		sprintf(infile3,"HaloRVList_halo.%.5d",nsteps[i]);
		npup= readdata(infile, infile1,infile2,infile3,&shup );
		debug(shup,npup, infile);
		npnow = npup;
		free(shnow); shnow = shup;
		MPI_Barrier(MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;

}

