#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<stddef.h>



int main (int argc, char **argv){
	int i,j,k;
	double size = 3150.L;
	FILE *fp, *wp;

	fp = fopen("fofhalo.z0.corrected.dat","r");
	wp = fopen("fofhalo.z0.dat","w");
	float mass,den,lden,hden,cden;
	int np;
	while((np=fscanf(fp,"%g %g %g %g %g \n",&mass,&den,&lden,&hden,&cden))>0){
		mass = mass/0.72;
		lden = lden/pow(0.72,3.L);
		hden = hden*pow(0.72,3.L);
		cden = cden*pow(0.72,3.L);
		fprintf(wp,"%g %g %g %g %g\n",mass,den,lden,hden,cden);
	}
	fclose(fp);
	fclose(wp);
}
