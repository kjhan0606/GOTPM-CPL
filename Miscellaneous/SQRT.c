#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>
struct timeval tv;
struct timezone tz;
float cputime0[100],cputime1[100];
#define TIMER_START(A) cputime0[A] = gettime();
#define TIMER_STOP(A)  cputime1[A] = gettime();
#define ELAPSED_TIME(A) (cputime1[A] - cputime0[A])


#define N 100000000

float gettime()
{
	static int startflag=1;
	static double tsecs0, tsecs1;
	
	if( startflag ) {
		(void) gettimeofday(&tv,&tz);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0e-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv,&tz);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;

	return ((float) (tsecs1-tsecs0));

}
double fsqrt(double number){
	long i;
	double x,y;
	const double f = 1.5L;
	if(number ==0.l) return 0.l;
	x = number*0.5L;
	y = number;
	i =*(long *)&y;
	i=0x5f3759df - (i>>1);
	y = *(double *) &i;
	y = y *(f-(x*y*y));
	y = y *(f-(x*y*y));
	return number*y;
}


float fsqrtf(float number){
	int i;
	float x,y;
	const float f = 1.5f;
	if(number ==0.f) return 0.f;
	x = number*0.5f;
	y = number;
	i =*(int *)&y;
	i=0x5f3759df - (i>>1);
	y = *(float *) &i;
	y = y *(f-(x*y*y));
	y = y *(f-(x*y*y));
	return number*y;
}

int main(){
	float a[N],b[N],c,d;
	int i,j,k;


	TIMER_START(1);
	for(i=0;i<N;i++){
		a[i] = 6./N*i;
		b[i] = fsqrtf(a[i]);
	}
	TIMER_STOP(1);
	TIMER_START(2);
	for(i=0;i<N;i++){
		a[i] = 6./N*i;
		b[i] = sqrtf(a[i]);
	}
	TIMER_STOP(2);
	printf("%g   <-----> %g\n",ELAPSED_TIME(1),ELAPSED_TIME(2));
	return 0;
}
