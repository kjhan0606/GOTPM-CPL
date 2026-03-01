#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  namu.exe 

#--- Objects to build

OBJS		= namu.o correl.o \
			pmigrate.o p_solver.o  pm_main_sub_pgi.o Memory2.o \
			timerutil.o  VarPM.o slicemigrate.o\
			savexzslice.o flagwriting.o force_spline.o\
			Treewalk.o  mydomaindecomposition.o   lightcone.o Pslicelightcone.o\
			kjhrw.o sub.lynx.camb.o ZeldovichInit.o dendump2disk.o\
			mysub.o header.o  determine_mpi_misc_param.o FindSeed.o \
			observerlightcone.obs0.o observerlightcone.obs1.o observerlightcone.obs2.o \
			observerlightcone.obs3.o observerlightcone.obs4.o observerlightcone.obs5.o \
			observerlightcone.obs6.o observerlightcone.obs7.o veltsc.o \
			indexing.o PencilBeamlightcone.o observerlightcone.obsCen.o\
			animate.o animate.obs0.o mkviewer.mod2.o ratint.o PreFoF.mod2.o fda4.o tsc.o \
 			kjhtpm.mod1.gpu.o  gpu.Treewalk.mod1.o Memory.gpu.o\
			# fda4_new.o fda4_split.o tsc_c.o

#			fmod.o fda4_new.o fda4_split.o tsc_c.o \
#indexing.o PMPI_Sendrecv_Replace.o  #kjhtpm.mod2.o
#hypersort.o\
#slicelightcone.final.o
#			kjhtpm.final.o indexing.o PMPI_Sendrecv_Replace.o  #kjhtpm.mod2.o
#			findindex.o maxextent.o glacial.o survey.o surveywrite2.o
#			memory.o	tpm.o domdecomp.o  kjhgetvisitors.o \
# 			tree.o  quickbin.o group.o groupgrav.o  \
#			kjhtpm.o Treewalk.mod6.o force_spline.mod2.o
#	hypersort.o \

OS = LINUX_GCC



#--- Targets

CAMBDIR = CAMB_Sources
#CAMBDIR = CAMB
SUBDIRS = $(CAMBDIR)

all: 
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) all);\
	done 
	$(MAKE) this

this: $(OBJS)
	rm -f $(PROGS) 
	$(MAKE) $(PROGS)
new: 
	$(MAKE) clean
	$(MAKE) all

clean: 
	rm -f $(OBJS)
	rm -f $(OBJPOW)
	rm -f $(OBJVIEW)
	rm -f $(PROGS)
	rm -f $(MKPOW)
	rm -f $(VIEWPOW)
	for subdir in $(SUBDIRS); do \
		(cd $$subdir && $(MAKE) clean);\
	done 

MKPOW = mkpower.exe
OBJPOW = mkpower.o
VIEWPOW = viewpower.exe
OBJVIEW = viewpower.o
MKPFILE = mkpfile.exe
MKTRACK = mktrack.exe
GPUTEST = gputest.exe

mkview:
	$(CC) -c -DNMEG=1000L mkviewer.mod2.c mktrack.c -g
	$(F90C) -o $(MKTRACK) mktrack.o mkviewer.mod2.o ratint.f Memory2.o -g
	

misc:
	rm -f $(VIEWPOW) 
	rm -f $(OBJVIEW) 
	rm -f $(MKPFILE) 
	rm -f $(MKTRACK) 
	$(F90C) -c viewpower.f90 -I./$(SUBDIRS) -O1  -mcmodel=medium
	$(F90C) -o $(VIEWPOW) viewpower.o -I./$(SUBDIRS) -L./$(CAMBDIR) -lcamb -O1 -mcmodel=medium
	$(F90C) -c mkpower.f90 -I./$(SUBDIRS) -O1  -mcmodel=medium
	$(CC) -c getsimparm.c -O1 -mcmodel=medium  -lm
	$(F90C) -o $(MKPOW) mkpower.o getsimparm.o header.o -I./$(SUBDIRS) -L./$(CAMBDIR) -lcamb -O1 -mcmodel=medium  
	$(CC) -o $(MKPFILE) mkpfile.c header.c -O1 -mcmodel=medium -lm
	$(MAKE) mkview

gputest:
	rm -rf {GPUTEST}
	gcc -c gputest.c   -g
	gcc -c timerutil.c   -g
	gcc -c Memory2.c   -g
	gcc -c force_spline.c   -g
	$(NVCC) -o $(GPUTEST) force_spline.o gputest.o Memory2.o gpu.Treewalk.mod1.cu \
		-lm -L$(CUDA)/lib64/ -L$(CUDA_SDK)/lib -I$(CUDA)/include \
		-I$(CUDA_SDK)/common/inc  timerutil.o  -g   -deviceemu

default: all

#$(MKPOW): $(OBJPOW)



include Rules.make
