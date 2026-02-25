#include<stdio.h>
#include<stdlib.h>
#include<stddef.h>
#include<string.h>
#include<math.h>
#include<mpi.h>
static int myid,nid;
int flagpsmeasure(float amax, float a, float astep,int step){
	FILE *fp;
	int istep;
	int saveflag;
	float red,redm,redp,redi;

	red = amax/a-1;
	redp = amax/(a-astep)-1;
	redm = amax/(a+astep)-1;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if(myid==0){
		if((fp=fopen("WriteSync&WholeDen.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if(fabs(red-redi) <= fabs(redp-redi) &&
						fabs(red-redi) <= fabs(redm-redi)){
					fclose(fp);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM_WORLD);
	if(step%20==1 || step<10 || fabs(amax-a) < 1.e-4) saveflag = 1;
	return saveflag;
}
int flagBinnedDen(float amax, float a, float astep){
	FILE *fp;
	int istep;
	int saveflag;
	float red,redm,redp,redi;

	red = amax/a-1;
	redp = amax/(a-astep)-1;
	redm = amax/(a+astep)-1;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if(myid==0){
		if((fp=fopen("BinnedDen.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if(fabs(red-redi) <= fabs(redp-redi) &&
						fabs(red-redi) <= fabs(redm-redi)){
					fclose(fp);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM_WORLD);

	return saveflag;
}
int flagPreFoF(float amax, float a, float astep, int stepcount){
	FILE *fp;
	int istep;
	int saveflag;
	float red,redm,redp,redi;

	red = amax/a-1;
	redp = amax/(a-astep)-1;
	redm = amax/(a+astep)-1;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if(myid==0){
		if((fp=fopen("PreFoF.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if(fabs(red-redi) <= fabs(redp-redi) &&
						fabs(red-redi) <= fabs(redm-redi)){
					fclose(fp);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM_WORLD);

	return saveflag;
}
int flagsyncpdata(float amax, float a, float astep){
	FILE *fp;
	int istep;
	int saveflag;
	float red,redm,redp,redi;

	red = amax/a-1;
	redp = amax/(a-astep)-1;
	redm = amax/(a+astep)-1;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if(myid==0){
		if((fp=fopen("WriteSyncPData.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if(fabs(red-redi) <= fabs(redp-redi) &&
						fabs(red-redi) <= fabs(redm-redi)){
					fclose(fp);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM_WORLD);
	return saveflag;
}
int flagwholeden(float amax, float a, float astep){
	FILE *fp;
	int istep;
	int saveflag;
	float red,redm,redp,redi;

	red = amax/a-1;
	redp = amax/(a-astep)-1;
	redm = amax/(a+astep)-1;
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	if(myid==0){
		if((fp=fopen("WriteWholeDen.flag","r"))){
			while(fscanf(fp,"%g",&redi)!= EOF){
				if(fabs(red-redi) <= fabs(redp-redi) &&
						fabs(red-redi) <= fabs(redm-redi)){
					fclose(fp);
					saveflag= 1;
					goto out;
				}
			}
			fclose(fp);
		}
	}
	saveflag = 0;
out: 
	MPI_Bcast(&saveflag,1,MPI_INT,0,MPI_COMM_WORLD);
	return saveflag;
}

int flagsuddenstop(int nowstep){
	int stopstep;
	int flag;
	FILE *fp;
	char flagcontinue[100];
	MPI_Comm_rank(MPI_COMM_WORLD,&myid);
	MPI_Comm_size(MPI_COMM_WORLD,&nid);

	flag= 0;
	if(myid==0){
		if((fp=fopen("Suddenstop.flag","r"))){
			fscanf(fp,"%d",&stopstep);
			if(stopstep == nowstep) {
				if(fscanf(fp,"%s",flagcontinue)!=EOF){
					if(strcmp(flagcontinue,"+")==0) flag = 2;
					else if(strcmp(flagcontinue,"-")==0) flag = 1;
					else flag = 3;
				}
				else flag= 1;
			}
			else flag= 0;
			fclose(fp);
		}
	}
	MPI_Bcast(&flag,1,MPI_INT,0,MPI_COMM_WORLD);
	return flag;
}

