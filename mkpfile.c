#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stddef.h>
#include<math.h>

#include "pmheader.h"
#include "params.h"


int main (int argc, char **argv){
	int i,j,k;
	FILE *wp;
	SimParameters sim;


	if(argc !=3){
		fprintf(stderr,"Please input as mkpfile.exe [outputfilename] [cosmology:WMAP3/WMAP5]\n");
		exit(99);
	}
	wp = fopen(argv[1],"w");
	{
		void mk_default_param(SimParameters *,char *);
		mk_default_param(&sim,argv[2]);
	}
	{
		void write_default_sim_parameter_file (FILE *, SimParameters);
		write_default_sim_parameter_file(wp,sim);
	}
	fclose(wp);


}
