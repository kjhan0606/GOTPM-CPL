#include<stdio.h>
#include<stdlib.h>

#include "prime.h"
#include "SimBox.h"

int getprimenumber(int number,PrimeNumber *prime){
	int i,j,k,denom,nprime;
	int res;


	nprime = 0;
	prime[nprime].factor = 0;
	while((number %2)==0){
		prime[nprime].prime = 2;
		prime[nprime].factor ++;
		number = number / 2;
	}
	if(prime[nprime].factor >0) nprime ++;
	prime[nprime].factor = 0;
	for(i=3;i<=number;i+=2){
		while((number % i) ==0){
			prime[nprime].prime = i;
			prime[nprime].factor ++;
			number = number/i;
		}
		if(prime[nprime].factor > 0) nprime ++;
		prime[nprime].factor = 0;
	}
	return nprime;
}

Range GetCorner(int np, basicparticletype *bp,
		int subgrpid, int nowgrpsize, int Ndivide, Range MotherBox,int axis){
	Range res;
	int nbin,i,j,k,ibin;
	float binsize,rbinsize;
	float *bin;
	/*
	res.min = MotherBox.min + (MotherBox.max-MotherBox.min)/Ndivide * subgrpid;
	res.max = MotherBox.min + (MotherBox.max-MotherBox.min)/Ndivide * (subgrpid+1);
	*/
	nbin = nowgrpsize *128;
	bin = (float*)Malloc(sizeof(float)*nbin,PPTR(bin));
	for(i=0;i<nbin;i++) bin[i] = 0;
	binsize = (MotherBox.max-MotherBox.min)/nbin;
	rbinsize = 1/binsize;
	for(i=0;i<np;i++){
		ibin = (bp[i].r[axis] - MotherBox.min)*rbinsize;
		bin[ibin] ++;
	}

	Free(bin);
	return res;
}

Box DivideBox(int nprime, PrimeNumber *prime, int NProc, int nowgrpid, Box box){
	int axis;
	int i,j,k,nowgrpsize;
	int subgrpsize,subgrpid;
	int nfamily,nneighbor;
	float lx,ly,lz;
	Box parentbox;

	nowgrpsize = NProc;
	for(i=nprime-1;i>=0;i--){
		while(prime[i].factor >0) {
			subgrpsize = nowgrpsize/prime[i].prime;
	
			lx = box.r[0].max - box.r[0].min;
			ly = box.r[1].max - box.r[1].min;
			lz = box.r[2].max - box.r[2].min;
			parentbox = box;
			subgrpid = nowgrpid/subgrpsize; /* To determine the local subgroup */
			if(lx >=ly && lx >= lz) axis = 0;
			else if(ly >= lx && ly >= lz) axis = 1;
			else axis = 2;
			box.r[axis] = GetCorner(subgrpid,nowgrpsize,prime[i].prime,box.r[axis],axis);
	
			nowgrpid = nowgrpid%subgrpsize;
			nowgrpsize = subgrpsize;
			prime[i].factor --;
		}
	}
	return box;
}



int main(int argc, char **argv){

	Box box;
	int nid,myid;
	int number,nprime,i;
	PrimeNumber *prime;

	printf("Please input number : ");
	fscanf(stdin,"%d",&nid);

	prime = (PrimeNumber*)malloc(sizeof(PrimeNumber)*(nid/2));

	nprime = getprimenumber(nid,prime);
	for(i=0;i<nprime;i++){
		printf("prime number is %d^%d\n",prime[i].prime,prime[i].factor);
	}


	for(myid=0;myid<nid;myid++){ 
		box.r[0].min = box.r[1].min = box.r[2].min = 0;
		box.r[0].max = box.r[1].max = box.r[2].max = 1024;
		nprime = getprimenumber(nid,prime);

		box = DivideBox(nprime,prime,nid,myid,box);
		printf("myid= %d has %g < x < %g && %g < y < %g && %g < z < %g\n",myid,
				box.r[0].min,box.r[0].max,box.r[1].min,box.r[1].max,box.r[2].min,box.r[2].max);
	}

	return 0;
}
