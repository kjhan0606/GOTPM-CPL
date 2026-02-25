#include <stdio.h>
#include <sys/time.h>
#include <sys/times.h>
#include <unistd.h>

struct timeval tv;
/*
struct timezone tz;
*/

#ifdef OLD
float gettime()
{
	float elapsedtime;
	struct tms rusage;
	
	times(&rusage);
	elapsedtime = (rusage.tms_utime + rusage.tms_stime)/((float) CLK_TCK);
	return(elapsedtime);
}
#else
float gettime()
{
	static int startflag=1;
	static double tsecs0, tsecs1;
	
	if( startflag ) {
		(void) gettimeofday(&tv,0);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0e-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv,0);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;

	return ((float) (tsecs1-tsecs0));

}
#endif

float wallclocktime()
{
	static int startflag=1;
	static double tsecs0, tsecs1, dt;

	if( startflag ) {
		(void) gettimeofday(&tv,0);
		tsecs0 = tv.tv_sec + tv.tv_usec*1.0e-6;
		startflag = 0;
	}
	(void) gettimeofday(&tv,0);
	tsecs1 = tv.tv_sec + tv.tv_usec*1.0e-6;
	dt = tsecs1-tsecs0;
	tsecs0 = tsecs1;
	if( dt <= 0 ) dt = 0;
	return (float) dt;
}
#include <time.h>
int wallclock()
{
	time_t tsecs0, dt;

	time(&tsecs0);
	dt = tsecs0;
	return ((int) dt);
}
