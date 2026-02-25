#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>

#include "pmheader.h"
#include "pm_common.h"
#include "flow.h"

void jwrite(char *, int );

void dumpmain(int iflagsyncpdata){
	if(halfstep.first==1 && halfstep.second==0) {
	/* Now positions and velocities of all simulation particles are 
	 * synchronized */
		if(iflagsyncpdata){
			char saveprefix[80],synchronizervprefix[80];
			strcpy(saveprefix,simpar.rvprefix);
			strcpy(synchronizervprefix,simpar.rvprefix);
			sprintf(simpar.rvprefix,"Sync%s",synchronizervprefix);
			jwrite(simpar.rvprefix,simpar.stepcount);
			strcpy(simpar.rvprefix,saveprefix);
		}
	}
}
