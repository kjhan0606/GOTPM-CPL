/* icc -o check check.c -DINDEX -DXYZDBL -g */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"
#include "mkmerger.h"



typedef struct MergeHistory{
	SavedHaloType a[75];
} MergeHistory;

int *nsteps;
float *reds;
float boxsize, hubble, npower, omep, omeplam, bias,astep,amax,anow,omepb;
int nx, nspace;
int msteps,*nsteps;


int main(int argc, char *argv[]){
	IDTYPE i,j,k;
	MergeHistory *bp;
	IDTYPE np;


	/*
	FILE *wp = fopen("MergingTree.dat","r");
	FILE *wp = fopen("sub2_MergingTree.dat","r");
	*/
	FILE *fp = fopen(argv[1],"r");
	FILE *wp = fopen(argv[2],"w");
	fread(&boxsize, sizeof(float), 1,fp);/* simulation box size in h^-1 Mpc */
	fread(&hubble, sizeof(float), 1,fp); /* hubble parameter, for example, h=0.72 */
	fread(&npower, sizeof(float), 1,fp); /* power spectral index */
	fread(&omep, sizeof(float), 1,fp); /* omega_matter at present */
	fread(&omepb, sizeof(float), 1,fp); /* Omega_baryon0 at present */
	fread(&omeplam, sizeof(float), 1,fp); /* Omega_lambda at present */
	fread(&bias, sizeof(float), 1,fp); /* bias factor or inverse of sigma_8 */
	fread(&astep, sizeof(float), 1,fp); /* simulation astep */
	fread(&msteps, sizeof(float), 1,fp); /* number of steps in the merging tree */
	nsteps = (int*)malloc(sizeof(int)*msteps);
	reds = (float*)malloc(sizeof(float)*msteps);
	fread(nsteps, sizeof(int), msteps,fp); /* step number save to merging tree */
	fread(reds, sizeof(float), msteps,fp); /* redshift info. save to merging tree */


	bp = (MergeHistory*)malloc(sizeof(MergeHistory)*1000000);
	size_t tcount,icount = 0;
	tcount = 0;
	int nblock = 1000000,iblock=0;
	float block[nblock];
	float rminmass,rmaxmass;
	int nbin = 30;
	double bin[nbin+1],massbin[nbin+1],err[nbin+1];
	rminmass = log10(2.7E11);
	rmaxmass = log10(5E15);
	float binstep = (rmaxmass-rminmass)/(float)nbin;
	for(i=0;i<=nbin;i++) {
		bin[i] = 0;
		massbin[i] = rminmass + (rmaxmass-rminmass)/(double) nbin * i;
	}
	long long gid= 0;
	long long ncount = 0;
	double boxsize = atof(argv[3]);

	while((np = fread(bp,sizeof(MergeHistory),1000000,fp))>0){
		for(i=0;i<np;i++){
			if(bp[i].a[msteps-1].mass > 0 && (gid == bp[i].a[msteps-1].majorglobalid)){
				float mass = log10(bp[i].a[msteps-1].fofhmass);
				int ntmp;
				ntmp = (mass-rminmass)/binstep;
				if(ntmp>=0 && ntmp < nbin) bin[ntmp] += 1.;
				ncount ++;
				
			}
			gid++;
		}
		printf("ruuning cumulative number is %ld\n",gid);fflush(stdout);
	}
	printf("Total FoF halos at z=0 is %ld\n",ncount);
	for(i=0;i<nbin;i++) {
                err[i] = sqrt(bin[i])/pow(boxsize,3);
                bin[i] = bin[i]/pow(boxsize,3);
        }

	for(i=0;i<nbin;i++){
		float binsize;
		float cmass, sigma,jenkins;
		binsize = pow(10,massbin[i+1]) - pow(10.L,massbin[i]);
		massbin[i] = 0.5*(pow(10.L,massbin[i]) + pow(10.L,massbin[i+1]));
		cmass = 0;
		for(j=i;j<nbin;j++) cmass += bin[j];
		printf("%lg %lg %lg %lg %lg\n",massbin[i],
                                bin[i]/binsize*massbin[i]*log(10.),
                                (bin[i]+err[i])/binsize*massbin[i]*log(10.),
                                (bin[i]-err[i])/binsize*massbin[i]*log(10.),
                                cmass);
                fprintf(wp,"%lg %lg %lg %lg %lg\n",massbin[i],
                                bin[i]/binsize*massbin[i]*log(10.),
                                (bin[i]+err[i])/binsize*massbin[i]*log(10.),
                                (bin[i]-err[i])/binsize*massbin[i]*log(10.),
                                cmass);

	}
	fclose(fp);
	fclose(wp);

}
