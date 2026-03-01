/* icc -o check check.c -DINDEX -DXYZDBL -g */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>
#include "merger.h"
#include "mkmerger.h"


#define Nsample 5


typedef struct MassHist{
	float minM,maxM;
	float red[75];
	double Mnow[75][1024];
	double psi[75][1024];
}MassHist;
MassHist msam[Nsample];



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


	msam[0].minM = 5e11;
	msam[0].maxM = 8e11;
	msam[1].minM = 1e12;
	msam[1].maxM = 3e12;
	msam[2].minM = 1e13;
	msam[2].maxM = 1e14;
	msam[3].minM = 1e14;
	msam[3].maxM = 5e14;
	msam[4].minM = 5e14;
	msam[4].maxM = 5e15;


	FILE *fp = fopen("MergingTree.dat","r");
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


	long long gid=0;


	for(i=0;i<Nsample;i++){
		for(j=0;j<msteps;j++){
			msam[i].red[j] = reds[j];
			for(k=0;k<1024;k++) {
				msam[i].Mnow[j][k] = 0;
				msam[i].psi[j][k] = 0;
			}
		}
	}

	float psimin,psimax;
	psimin = 0;
	psimax = 5;


	bp = (MergeHistory*)malloc(sizeof(MergeHistory)*100000);
	long long icount = 0;

	while((np = fread(bp,sizeof(MergeHistory),100000,fp))>0)
	{
		for(i=0;i<np;i++){
			if(gid != bp[i].a[msteps-1].majorglobalid || bp[i].a[msteps-1].mbp <=0){
				gid ++;
				continue;
			}
			float fmass = bp[i].a[msteps-1].mass;
			for(j=0;j<Nsample;j++){
				if(fmass >= msam[j].minM && fmass < msam[j].maxM) {
					break;
				}
			}
			if(j == Nsample) {
				gid ++;
				continue;
			}
			int jsample = j;
			icount ++;
			for(j=1;j<msteps;j++){
				if(bp[i].a[j].mbp >0){
					/* if bp[i].a[j].mbp ==0 --> mbp info is null at this redshift. Perhaps not yet formed */
					k = (bp[i].a[j].mass/fmass)*1024;
					int kk = (-log(bp[i].a[j].mass/fmass)/reds[j] -psimin)/(psimax-psimin) * 1024;
					if(k<1024) msam[jsample].Mnow[j][k] += 1;
					if(kk<1024 && kk >=0) msam[jsample].psi[j][kk] += 1;
				}
			}
			gid ++;
		}
	}
	fclose(fp);

	printf("Total %ld halo counted\n",icount);


	FILE *wp = fopen("mergingdist.out","w");
	for(i=0;i<Nsample;i++){
		fprintf(wp,"\n\n # sample %g <= M < %g\n\n",msam[i].minM, msam[i].maxM);
		for(j=0;j<msteps;j++){
			long isum=0;
			for(k=0;k<1024;k++) {
				isum += msam[i].Mnow[j][k];
			}
			if(isum>=10) {
				long zsigma = isum/2;
				long psigma = 0.683*isum;
				long msigma = 0.317*isum;
				float zsig,psig,msig;
				long lsum = 0;
				int mflag,pflag,cflag;
				mflag=pflag=cflag = 0;
				for(k=0;k<1024;k++){
					lsum += msam[i].Mnow[j][k];
					if(mflag ==0 && lsum > msigma) {
						msig = k/1024.;
						mflag = 1;
					} 
					if(cflag ==0 && lsum > zsigma) {
						zsig = k/1024.;
						cflag = 1;
					} 
					if(pflag ==0 && lsum > psigma) {
						psig = k/1024.;
						pflag = 1;
					} 
				}
				fprintf(wp,"%g %g %g %g ",msam[i].red[j],msig,zsig,psig);

				lsum = 0;
				mflag=pflag=cflag = 0;
				for(k=0;k<1024;k++){
					lsum += msam[i].psi[j][k];
					if(mflag ==0 && lsum > msigma) {
						msig = k/1024. *(psimax-psimin) + psimin;
						mflag = 1;
					} 
					if(cflag ==0 && lsum > zsigma) {
						zsig = k/1024. *(psimax-psimin) + psimin;
						cflag = 1;
					} 
					if(pflag ==0 && lsum > psigma) {
						psig = k/1024. *(psimax-psimin) + psimin;
						pflag = 1;
					} 
				}
				fprintf(wp," %g %g %g\n",msig,zsig,psig);
			}
		}
	}
	fclose(wp);
}
