#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#define IMOD(A,B) ((A) - ((A)/(B))*(B))

int main(){
	int i,j,k;

	for(i=0;i<8;i++){
		j = IMOD(i,2);
		k = (IMOD(i,4)/2);
		printf("%d %d %d\n",i,j,k);
	}
	for(i=0;i<8;i++){
		j = (i%2);
		k = ((i%4)/2);
		printf("%d %d %d\n",i,j,k);
	}
}
