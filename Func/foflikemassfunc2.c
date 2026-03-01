#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "link.h"
#define PI 3.141592655L
#define MAX_HALO_NUM 10000000
#define EPS 0.05

/*
#define MIN(A,B) ((A) > (B) ? (B):(A))
#define MAX(A,B) ((A) < (B) ? (B):(A))
*/
typedef struct fofhalo{
	float mass;
} fofhalo;
fofhalo *mhalo;
float boxsize;
int nummhalo;
float R1[3][3],R2[3][3],D1[3][3],D2[3][3];
float T1[3],T2[3];
float A[3][3],B[3][3];
#define N 20
#define M 40
#define NM 800
void SEMI_FoF_LINK(Halo *halo,int nummass){
	float tmpmass;
	int phalo[10000],nphalo;
	int i,j,k;
	float tmpdist;
	int ntmp;
	float xtmp,ytmp,ztmp;
	float dist;
	float xp[3][NM];
	float px[3][NM];
	float pxtmp[3][NM],xptmp[3][NM];
	mhalo = (fofhalo *) malloc(sizeof(fofhalo)*nummass);
	nummhalo = 0;
	ntmp = 0;
	for(i=0;i<N;i++)
		for(j=0;j<M;j++){
			xp[0][ntmp] = sin(PI/N * (float)i)*cos(2.*PI/M*(float)j);
			xp[1][ntmp] = sin(PI/N * (float)i)*sin(2.*PI/M*(float)j);
			xp[2][ntmp] = cos(PI/N * (float)i);
			ntmp++;
		}
	for(j=0;j<3;j++){
		D1[0][j] = 0;
		D1[1][j] = 0;
		D1[2][j] = 0;
		D2[0][j] = 0;
		D2[1][j] = 0;
		D2[2][j] = 0;
	}

	for(i=0;i<nummass;i++){
		if(halo[i].mpeak == 1){
			halo[i].included = 0;
			mhalo[nummhalo].mass = halo[i].mass;
			nummhalo ++;
		}
		else if(halo[i].included == -1 && halo[i].mpeak > 1 ){
			if(fmod(i,10) == 0) printf("%d complete in %d\n",i,nummass);
			phalo[0] = i;
			halo[i].included = 0;
			nphalo = 1;
			for(k=0;k<nphalo;k++){
				int ii;
				ii = phalo[k];
				for(j=0;j<3;j++){
					R1[0][j] = halo[ii].a[j];
					R1[1][j] = halo[ii].b[j];
					R1[2][j] = halo[ii].c[j];
				}
				D1[0][0] = halo[ii].r;
				D1[1][1] = halo[ii].r*halo[ii].q;
				D1[2][2] = halo[ii].r*halo[ii].s;
				/*
				D1[0][0] = halo[ii].r;
				D1[1][1] = halo[ii].r;
				D1[2][2] = halo[ii].r;
				*/
				T1[0] = halo[ii].x;
				T1[1] = halo[ii].y;
				T1[2] = halo[ii].z;
				{
					int ii,jj,kk;
					for(ii=0;ii<3;ii++)
						for(jj=0;jj<3;jj++){
							B[ii][jj] = 0;
							for(kk=0;kk<3;kk++){
								B[ii][jj] += R1[ii][kk]*D1[kk][jj];
							}
						}
					for(ii=0;ii<N*M;ii++){
						for(jj=0;jj<3;jj++){
							px[jj][ii] = 0;
							for(kk=0;kk<3;kk++){
								px[jj][ii] += B[jj][kk]*xp[kk][ii];
							}
						}
					}
					for(ii=0;ii<N*M;ii++){
						px[0][ii] = px[0][ii] + T1[0];
						px[1][ii] = px[1][ii] + T1[1];
						px[2][ii] = px[2][ii] + T1[2];
					}
				}
	
				for(j=i+1;j<nummass;j++){
					if(halo[j].indx != halo[i].indx) break;
					if(halo[j].included == -1){
						float xxtmp,yytmp,zztmp;
						xtmp = halo[j].x-halo[phalo[k]].x;
						ytmp = halo[j].y-halo[phalo[k]].y;
						ztmp = halo[j].z-halo[phalo[k]].z;
						T2[0] = halo[j].x;
						T2[1] = halo[j].y;
						T2[2] = halo[j].z;
						if(xtmp > boxsize*0.5) {
							xtmp = boxsize - xtmp;
							T2[0] -= boxsize;
						}
						else if(xtmp < -boxsize*0.5){
							xtmp = boxsize + xtmp;
							T2[0] += boxsize;
						}
						if(ytmp > boxsize*0.5) {
							ytmp = boxsize - ytmp;
							T2[1] -= boxsize;
						}
						else if(ytmp < -boxsize*0.5){
							ytmp = boxsize + ytmp;
							T2[1] += boxsize;
						}
						if(ztmp > boxsize*0.5) {
							ztmp = boxsize - ztmp;
							T2[2] -= boxsize;
						}
						else if(ztmp < -boxsize*0.5){
							ztmp = boxsize + ztmp;
							T2[2] += boxsize;
						}
						tmpdist = xtmp*xtmp+ytmp*ytmp+ztmp*ztmp;
						tmpdist = sqrt(tmpdist);
						if(tmpdist <= halo[phalo[k]].r+halo[j].r){
							int ii,jj,kk;
							for(ii=0;ii<3;ii++){
								R2[0][ii] = halo[j].a[ii];
								R2[1][ii] = halo[j].b[ii];
								R2[2][ii] = halo[j].c[ii];
							}
							D2[0][0] = 1./halo[j].r;
							D2[1][1] = 1./halo[j].r/halo[j].q;
							D2[2][2] = 1./halo[j].r/halo[j].s;
							/*
							D2[0][0] = 1./halo[j].r;
							D2[1][1] = 1./halo[j].r;
							D2[2][2] = 1./halo[j].r;
							*/
							for(ii=0;ii<3;ii++)
								for(jj=0;jj<3;jj++){
									A[ii][jj] = 0;
									for(kk=0;kk<3;kk++) {
										A[ii][jj] += D2[ii][kk]*R2[kk][jj];
									}
								}
							for(ii=0;ii<NM;ii++){
								for(jj=0;jj<3;jj++){
									pxtmp[jj][ii] = px[jj][ii]-T2[jj];
								}
							}
							for(ii=0;ii<NM;ii++){
								for(jj=0;jj<3;jj++){
									xptmp[jj][ii] = 0;
									for(kk=0;kk<3;kk++){
										xptmp[jj][ii] += 
											A[jj][kk] * pxtmp[kk][ii];
									}
								}
							}
							for(ii=0;ii<NM;ii++){
								xtmp = xptmp[0][ii];
								ytmp = xptmp[1][ii];
								ztmp = xptmp[2][ii];
								tmpdist = sqrt(xtmp*xtmp+ytmp*ytmp+ztmp*ztmp);
								if(tmpdist <= 1.-EPS) break;
							}
							if(ii < NM){
								halo[j].included = 0;
								phalo[nphalo] = j;
								nphalo++;
							}
						}
					}
				}
			}
			tmpmass = 0;
			for(k=0;k<nphalo;k++){
				tmpmass += halo[phalo[k]].mass;
			}
			mhalo[nummhalo].mass = tmpmass;
			nummhalo ++;
		}
	}
};

int main(int argc, char **argv){
	int i,j,k;
	char filename[80];
	float mass[MAX_HALO_NUM];
	int nummass;
	float minmass=1.E23,maxmass=-1.E23;
	float *massbin,*numbin, *err;
	float binstep;
	Halo *halo;
	FILE *fp,*wp;
	char *outfile;
	int nbin;
	float hubble;
	float x,y,z,index,mpeak,peakindex;

	if(argc != 6) {
		fprintf(stderr,
				"input : massfunc2 [filename] [# of bin] [boxsize] [h] [outfile]\n");
		exit(0);
	}
	fp = fopen(argv[1],"r");
	nbin = atoi(argv[2]);
	boxsize = atof(argv[3]);
	hubble = atof(argv[4]);
	outfile = argv[5];
	halo = (Halo *)malloc(sizeof(Halo)*MAX_HALO_NUM);
	i = 0;
	/*
	while(fscanf(fp,"%f %f %f %f %d %d %d\n",mass+i,&x,&y,&z,
				&index,&mpeak,&peakindex) != EOF){
				*/
	while(fscanf(fp,"%f %f %f %f %d %d %d %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g %g\n",
			&(halo[i].mass),&(halo[i].x),&(halo[i].y),&(halo[i].z),
		    &(halo[i].indx),&(halo[i].mpeak),&(halo[i].peakindx),
		    &(halo[i].spin),&(halo[i].pv),&(halo[i].ang),
		    &(halo[i].r),&(halo[i].q),&(halo[i].s),
		    &(halo[i].a[0]),&(halo[i].a[1]),&(halo[i].a[2]),
		    &(halo[i].b[0]),&(halo[i].b[1]),&(halo[i].b[2]),
		    &(halo[i].c[0]),&(halo[i].c[1]),&(halo[i].c[2]),
		    &(halo[i].rot[0]),&(halo[i].rot[1]),&(halo[i].rot[2]))!= EOF){

		/* multiply by hubble  because massfunc1.c outputs 
		 * mass including hubble constant */
		halo[i].mass *= hubble;
		minmass = MIN(minmass,halo[i].mass);
		maxmass = MAX(maxmass,halo[i].mass);
		halo[i].included = -1;
		/*
		*(mass+i) = *(mass+i)*hubble;
		minmass = MIN(minmass,*(mass+i));
		maxmass = MAX(maxmass,*(mass+i));
		*/
		i++;
	}
	nummass = i;
	fclose(fp);
	/*
	 * Merging halos within each virial radius(circular implementation)
	 * to higher mass halo
	 */
	SEMI_FoF_LINK(halo,nummass);

	printf("minum %g and maximum %g\n",minmass,maxmass);
	massbin = (float *) malloc(sizeof(float)*(nbin+1));
	numbin = (float *) malloc(sizeof(float)*(nbin+1));
	err = (float *) malloc(sizeof(float)*(nbin+1));
	minmass = log10(minmass);
	maxmass = log10(maxmass);
	for(i=0;i<=nbin;i++){
		massbin[i] = minmass + (maxmass-minmass)/(float) nbin *i;
		numbin[i] = 0;
	}
	binstep = (maxmass-minmass)/(float) nbin;
	for(i=0;i<nummhalo;i++){
		float tmpmass;
		int ntmp;
		tmpmass = log10(mhalo[i].mass);
		ntmp = (tmpmass-minmass)/binstep;
		numbin[ntmp] += 1;
	}
	/* centering of bin's to take account of the [dex] */
	wp = fopen(outfile,"w");
	for(i=0;i<nbin;i++) {
		err[i] = sqrt(numbin[i])/pow(boxsize,3);
		numbin[i] = numbin[i]/pow(boxsize,3);
	}
	for(i=0;i<nbin;i++){
		float binsize;
		float cmass;
		float sigma,jenkins;
		binsize = pow(10,massbin[i+1]) - pow(10,massbin[i]);
		massbin[i] = (pow(10,massbin[i]) + pow(10,massbin[i+1]))*0.5;
		cmass = 0;
		for(j=i;j<nbin;j++){
			cmass += numbin[j];
		}
		printf("%g %g %g %g %g\n",massbin[i],
				numbin[i]/binsize*massbin[i]*log(10.),
				(numbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(numbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
		fprintf(wp,"%g %g %g %g %g\n",massbin[i],
				numbin[i]/binsize*massbin[i]*log(10.),
				(numbin[i]+err[i])/binsize*massbin[i]*log(10.),
				(numbin[i]-err[i])/binsize*massbin[i]*log(10.),
				cmass);
	}
	fclose(wp);
	return 0;
}
