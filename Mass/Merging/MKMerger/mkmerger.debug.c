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

IDTYPE globalid = 0;


#define NSTEPS 75
int nsteps[NSTEPS] ={136, 165, 183, 205, 234, 250, 270, 292, 318, 
	321, 353, 385, 415, 444, 473, 501,  555, 632, 657, 681, 
	706, 730, 754, 801, 824, 916, 938, 961, 983, 1006, 1028, 1051, 
	1073, 1095, 1118, 1140, 1162, 1184, 1207, 1229, 1244, 1252, 1274, 1297, 1319, 
	1328, 1342, 1365, 1388, 1410, 1434, 1457, 1480, 1503, 1527, 1535, 1550, 
	1574, 1598, 1622, 1646, 1670, 1694, 1738, 1744, 1769, 1794, 1819, 1844, 
	1870, 1896, 1922, 1948, 1974, 2001};


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
	
	printf("P%d has np = %ld : shalo = %p\n",myid,np,shalo);
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
		shalo[i].basichalo.nowhid = (shalo[i].nowhid = (shalo[i].globalid = -1));
		shalo[i].basichalo.mbp = bp[i].mbp;
		shalo[i].basichalo.nowhid = shalo[i].nowhid = bp[i].nowhid;
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
void dumpdata(SimplifiedHalo *shnow, IDTYPE npnow,char *outfile, int istep){
	IDTYPE i,j,k,nshift,gidsize;
	FILE *wp;


	int isend, iget ; 
	isend = iget =1; 
	int WGroupSize = WGROUPSIZE; 
	int src = myid-1; 
	int tgt = myid+1; 
	int itag = 10;
	MPI_Status status;
	if(myid==0) printf("\n\n\n #################################\n\n");
#ifdef PlEASE_CHECK
	/* If the multiple open-to-write mode can be applied */
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);

#else
	for(j=0;j<nid;j++)
	{
		MPI_Barrier(MPI_COMM_WORLD);
		if(myid==j)
#endif
		{
			if(npnow >0) {
				wp = fopen(outfile,"rw+");
				gidsize = sizeof(SavedHaloType)*NSTEPS;
			
				IDTYPE sgid = shnow[0].globalid;

#ifdef SLOW
				fseek(wp,initshift + sgid * gidsize +istep*sizeof(SavedHaloType), SEEK_SET);
				for(i=0;i<npnow;i++){
					nshift  = gidsize * (shnow[i].globalid - sgid);
					fseek(wp,nshift,SEEK_CUR);
					fwrite(&(shnow[i].basichalo),sizeof(SavedHaloType),1,wp);
					sgid = shnow[i].globalid;
				}
#else
				long offset ;
				if((offset=fseek(wp,initshift + sgid * gidsize , SEEK_SET))){
					fprintf(stderr,"P%d Error in shifting the file pointer %ld with initshift= %ld\n",myid, offset,
							initshift);
				}
				else {
					if(ftell(wp) != initshift + sgid*gidsize){
						fprintf(stderr,"P%d has problem in offset %ld ::: %ld\n",myid, initshift+sgid*gidsize, ftell(wp));
					}
				}
				IDTYPE nbites = sizeof(SavedHaloType)*(shnow[npnow-1].globalid - sgid + 1);
				IDTYPE nblockid = 400000;
				IDTYPE blocksize = sizeof(SavedHaloType)*nblockid;
				for(i=0;i<npnow;i+= nblockid){
					IDTYPE mp = MIN((npnow-i), nblockid);
					SavedHaloType *buff = (SavedHaloType *)calloc(mp*NSTEPS, sizeof(SavedHaloType));
					size_t mmp = fread(buff, sizeof(SavedHaloType)*NSTEPS, mp, wp);
					for(k=0;k<npnow;k++){
						long long ioffset = shnow[k].globalid - (sgid + i);
						if(ioffset <0) continue;
						else if(ioffset >= mp) break;
						else {
							SavedHaloType *tmp = buff + (ioffset*NSTEPS + istep);
							*tmp = shnow[k].basichalo;
						}
	
					}
					fseek(wp,-mmp*gidsize,SEEK_CUR);
					fwrite(buff, sizeof(SavedHaloType)*NSTEPS, mp,wp);
					free(buff);
				}
#endif
				fclose(wp);
			}
			printf("P%d has written data %ld\n",myid,npnow);
		}
#ifndef PLEASE_CHECK
	}
#else
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
#endif
}
void advdumpdata(SimplifiedHalo *shnow, IDTYPE npnow,char *outfile){
	IDTYPE i,j,k,nshift,gidsize;
	FILE *wp;


	int isend, iget ; 
	isend = iget =1; 
	int WGroupSize = WGROUPSIZE; 
	int src = myid-1; 
	int tgt = myid+1; 
	int itag = 10;
	MPI_Status status;
	if(myid==0) {
		printf("\n\n\n #################################\n\n");
		wp = fopen(outfile,"w");
		fclose(wp);
	}
#ifdef PLEASE_CHECK
	/* If the multiple open-to-write mode can be applied */
	if(RANKINGROUP(myid,WGroupSize) != 0 ) 
		MPI_Recv(&iget,1,MPI_INT,src,itag,MPI_COMM_WORLD,&status);

#else
	for(j=0;j<nid;j++)
	{
		if(myid==j)
#endif
		{
			if(myid==0) wp = fopen(outfile,"w");
			else wp = fopen(outfile,"rwa+");
			if(npnow >0) {
				gidsize = sizeof(SavedHaloType);
			
				IDTYPE sgid = shnow[0].globalid;

#ifdef SLOW
				fseek(wp,initshift + sgid * gidsize , SEEK_SET);
				for(i=0;i<npnow;i++){
					nshift  = gidsize * (shnow[i].globalid - sgid);
					fseek(wp,nshift,SEEK_CUR);
					fwrite(&(shnow[i].basichalo),sizeof(SavedHaloType),1,wp);
					sgid = shnow[i].globalid;
				}
#else
				long offset ;
				if((offset=fseek(wp,0 + sgid * gidsize , SEEK_SET))){
					fprintf(stderr,"P%d Error in shifting with %ld the file pointer %ld\n",myid, sgid*gidsize,offset);
					exit(99);
				}
				else {
					if(ftell(wp) != 0 + sgid*gidsize){
						fprintf(stderr,"P%d has problem in offset %ld ::: %ld\n",myid, 0+sgid*gidsize, ftell(wp));
						exit(99);
					}
				}
				IDTYPE nblockid = 400000;
				IDTYPE blocksize = sizeof(SavedHaloType)*nblockid;
				SavedHaloType *buff = (SavedHaloType *)calloc(nblockid, sizeof(SavedHaloType));
				for(k=0;k<nblockid;k++){
					buff[k].mbp = buff[k].nowhid = buff[k].majorglobalid =  -1;
				}
				IDTYPE jj;
				for(i=0;i<npnow;i++){
					jj = shnow[i].globalid - sgid;
					if(jj >= nblockid) {
						fwrite(buff, sizeof(SavedHaloType), nblockid,wp);
						for(k=0;k<nblockid;k++){
							buff[k].mbp = buff[k].nowhid = buff[k].majorglobalid =  -1;
						}
						sgid = shnow[i-1].globalid + 1;
					}
					jj = shnow[i].globalid - sgid;
					buff[jj] = shnow[i].basichalo;
				}
				fwrite(buff, sizeof(SavedHaloType), jj+1,wp);
				free(buff);
/*
				IDTYPE runninggid = sgid;
				for(i=0;i<npnow;i+= nblockid){
					IDTYPE mp = MIN((npnow-i), nblockid);
					IDTYPE mmp = shnow[mp-1].globalid - shnow[i].globalid + 1;
					SavedHaloType *buff = (SavedHaloType *)calloc(mmp, sizeof(SavedHaloType));
					for(k=i;k<npnow;k++){
						long long ioffset = shnow[k].globalid - (shonow[i].globalid);
						if(ioffset <0) continue;
						else if(ioffset >= mmp) break;
						else {
							buff[ioffset] = shnow[k].basichalo;
						}
	
					}
					if(k <npnow) runninggid = shnow[k].globalid;
					fwrite(buff, sizeof(SavedHaloType), mmp,wp);
					free(buff);
				}
*/
#endif
			}
			fclose(wp);
			printf("P%d has written data %ld\n",myid,npnow);
		}
		MPI_Barrier(MPI_COMM_WORLD);
#ifndef PLEASE_CHECK
	}
#else
	if(GROUPID(myid,WGroupSize) == GROUPID(tgt,WGroupSize) && tgt < nid) 
		MPI_Send(&isend,1,MPI_INT,tgt,itag,MPI_COMM_WORLD);
#endif
}


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


void dumpheader(char *outfile){
	FILE *wp;

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
}
void debug(SimplifiedHalo *arr, IDTYPE np){
	return;
	IDTYPE i;
	for(i=0;i<np;i++){
		if(arr[i].nowaid == 0 && arr[i].upaid==0){
			printf("P%d errror hererere %ld\n",myid,i);
		}
	}
}
IDTYPE debug1(SimplifiedHalo *arr, IDTYPE np){
	IDTYPE i;
	IDTYPE ncount = 0;
	for(i=0;i<np;i++){
		if(arr[i].globalid >= 0){
			ncount ++;
		}
	}
	return ncount;
}

int main(int argc, char **argv){
	IDTYPE i,j,k;
	IDTYPE np1,np2,ip;
	TrHalo *nmbp,imbp;
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

	/* Setting initial global id to the first slice data */
	if(myid!=0) MPI_Recv(&globalid,sizeof(IDTYPE),MPI_BYTE,myid-1,0,MPI_COMM_WORLD,&status);
	for(i=0;i<npnow;i++){
		shnow[i].globalid = globalid;
		globalid++;
	}
	if(myid != nid-1) MPI_Send(&globalid,sizeof(IDTYPE),MPI_BYTE,myid+1,0,MPI_COMM_WORLD);
	MPI_Bcast(&globalid,sizeof(IDTYPE),MPI_BYTE,nid-1,MPI_COMM_WORLD);


	determinemajorglobalid(shnow, npnow);
	if(myid==0) dumpheader(outfile);



	initshift = sizeof(float)*8 + sizeof(int)*1 + (sizeof(int)+sizeof(float))*NSTEPS;





#ifdef NONADVSAVE
	dumpdata(shnow, npnow, outfile,0);
#else
	{
		char out1[190]; sprintf(out1,"%s.00000",outfile);
		advdumpdata(shnow, npnow, out1);
	}
#endif

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
		debug(shup,npup);
		/*
		IDTYPE nnexthid = readhidmass(infile2,&hidmass);
		*/

		/* assign globalid to the next-step mbp data */
		{
			IDTYPE min,max,gmin,gmax;
			getminmaxofupnowaid(shnow,npnow, shup,npup, &gmin, &gmax);
			min = gmin + myid*((gmax-gmin-1)/nid+1);
			max = gmin + (myid+1)*((gmax-gmin-1)/nid+1);
			min = MIN(min,gmax);
			max = MIN(max,gmax);
			if(max == gmax) max = gmax + 1;
			MyFixedRangePqsort(shnow,npnow,SimplifiedHalo,IDTYPE,upaid,min,max,
					MPI_COMM_WORLD);
			MyFixedRangePqsort(shup,npup,SimplifiedHalo,IDTYPE,nowaid,min,max,
					MPI_COMM_WORLD);
			k = 0;
			IDTYPE icount = 0;
			for(j=0;j<npnow;j++){
				for(;k<npup;k++){
					if(shnow[j].upaid == shup[k].nowaid){
						shup[k].globalid = shnow[j].globalid;
						if(shnow[j].basichalo.mbp != shup[k].basichalo.mbp){
							icount ++;
							printf("P%d get strange %ld : %ld (%ld) ||||| %ld %ld : %ld %ld\n",myid, j,k,npnow,
									shnow[j].basichalo.mbp,shup[k].basichalo.mbp,
									shnow[j].nowaid, shup[k].nowaid);
						}
					}
					else if(shup[k].nowaid > shnow[j].upaid) break;
				}
			}

			if(icount) {
				printf("P%d has icount = %ld\n",myid,icount);
				exit(99);
			}
			IDTYPE nlink = debug1(shup,npup);
			{
				float frac = (float)nlink / (float) npnow;
				printf("P%d has linking success %ld : %ld :::; %g\n",myid,nlink,npnow,frac);
			}
			IDTYPE tnlink,tnpnow;
			MPI_Reduce(&nlink,&tnlink,1,MPI_IDTYPE, MPI_SUM,0,MPI_COMM_WORLD);
			MPI_Reduce(&npnow,&tnpnow,1,MPI_IDTYPE, MPI_SUM,0,MPI_COMM_WORLD);
			if(myid==0){	
				float frac = (float)tnlink / (float) tnpnow;
				printf("Total linking success is %ld : %ld :::; %g\n",tnlink,tnpnow,frac);
			}
			MPI_Barrier(MPI_COMM_WORLD);


		}

		{
			IDTYPE min,max,tmin,tmax;
			min = INT_MAX;
			max = INT_MIN;
			for(j=0;j<npup;j++){
				if(shup[j].nowhid > max) max = shup[j].nowhid;
				if(shup[j].nowhid < min) min = shup[j].nowhid;
			}
			long lmin,ltmin,lmax,ltmax;
			lmin = min;
			MPI_Reduce(&lmin,&ltmin,1,MPI_LONG,MPI_MIN,0,MPI_COMM_WORLD);
			lmax = max;
			MPI_Reduce(&lmax,&ltmax,1,MPI_LONG,MPI_MAX,0,MPI_COMM_WORLD);
			if(myid==0) {
				min = ltmin;
				max = ltmax;
			}
			MPI_Bcast(&min, 1, MPI_IDTYPE, 0, MPI_COMM_WORLD);
			MPI_Bcast(&max, 1, MPI_IDTYPE, 0, MPI_COMM_WORLD);
			IDTYPE gmin,gmax;
			gmin = min; gmax = max;
			min = gmin + myid*((gmax-gmin-1)/nid+1);
			max = gmin + (myid+1)*((gmax-gmin-1)/nid+1);
			min = MIN(min,gmax);
			max = MIN(max,gmax);
			if(myid==nid-1) max = gmax + 1;
			MyFixedRangePqsort(shup,npup,SimplifiedHalo,IDTYPE,nowhid,min,max,
					MPI_COMM_WORLD);
			determinemajorglobalid(shup, npup);
		}

		debug(shup,npup);



		/* assign global id which is the array index that will be save to a file */
		if(myid!=0) MPI_Recv(&globalid,1,MPI_IDTYPE,myid-1,0,MPI_COMM_WORLD,&status);
		for(k=0;k<npup;k++){
			if(shup[k].globalid<0){
				shup[k].globalid = globalid;
				globalid++;
			}
		}
		if(myid != nid-1) MPI_Send(&globalid,1,MPI_IDTYPE,myid+1,0,MPI_COMM_WORLD);
		MPI_Bcast(&globalid,1,MPI_IDTYPE,nid-1,MPI_COMM_WORLD);

		/* sort data in the order of globalid */
		{
			IDTYPE gmin,gmax,min,max;
			gmin = 0;
			gmax = globalid+1;
			min = gmin + myid*((gmax-gmin-1)/nid+1);
			max = gmin + (myid+1)*((gmax-gmin-1)/nid+1);
			min = MIN(min,gmax);
			max = MIN(max,gmax);
			if(min > gmax) min = gmax;
			if(max > gmax) max = gmax;
			/*
			MyConstRangePqsort(shup,npup,sizeof(SimplifiedHalo),SimplifiedHalo,IDTYPE,globalid,min,max,
					MPI_COMM_WORLD);
					*/
			MyFixedRangePqsort(shup,npup,SimplifiedHalo,IDTYPE,globalid,min,max,
					MPI_COMM_WORLD);
		}
		debug(shup,npup);

#ifdef NONADVSAVE
		dumpdata(shup, npup, outfile,i);
#else
		{
			char out1[190]; sprintf(out1,"%s.%.5d",outfile,(int)i);
			advdumpdata(shup, npup, out1);
		}
#endif
		debug(shup,npup);

		npnow = npup;
		free(shnow); shnow = shup;
	}

#ifdef NONADVSAVE
#else
/*
	if(myid==0) MergeFiles(outfile,globalid);
*/
#endif

	MPI_Finalize();
	return 0;

}

