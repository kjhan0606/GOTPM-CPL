#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  namu.exe 


include sRules.make


#--- Objects to build

OBJS		= main.o namu.zeldovich.o namu.2lpt.o addTwoOrderLPT.o TwoLPT.o TwoOrderLPT.o\
			pmmain.o dumpmain.o treemain.o \
			correl.o \
			p_solver.o  Memory2.o \
			timerutil.o  VarPM.o slicemigrate.o\
			savexzslice.o flagwriting.o force_spline.o\
			Treewalk.o  \
			kjhrw.o sub.lynx.camb.o addZeldovichInit.o ZeldovichInit.o dendump2disk.o\
			$(MYSUB).o header.o  determine_mpi_misc_param.o FindSeed.o \
			indexing.o  pcorr.o \
			ratint.o fda4.o tsc.o FLRW.o NR.o\
			fft.o pm_main_sub_fftw3.o veltsc.o\
			tsc_omp2.o pmigrate.kisti.omp.o mydomaindecomposition.omp.o kjhtpm.mod1.gpu.omp.o PreFoF.mod2.omp.o

#			Disklightcone.o Bosslightcone.o\
			lightcone.o Pslicelightcone.o\
 			observermain.mod2.o  Gott.mod1.o \
			animate.o animate.obs0.o mkviewer.mod2.o
	#
#
#
#
#pmigrate.kisti.o mydomaindecomposition.o kjhtpm.mod1.gpu.o PreFoF.mod2.o
#


		#	pm_main_sub_pgi_c.o  \
# gpu.Treewalk.mod2.o Memory.gpu.o\
#			observermain.o \
#			observerlightcone.obs0.o observerlightcone.obs1.o observerlightcone.obs2.o \
#			observerlightcone.obs3.o observerlightcone.obs4.o observerlightcone.obs5.o \
#			observerlightcone.obs6.o observerlightcone.obs7.o\
#			PencilBeamlightcone.o observerlightcone.obsCen.o\
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

#CAMBDIR = CAMB_Sources
#CAMBDIR = CAMB
#SUBDIRS = $(CAMBDIR)

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
PCORRTEST = pcorrtest.exe

mkview:
	$(CC) -c -DNMEG=1000L mkviewer.mod2.c mktrack.c -g
	$(F90C) -o $(MKTRACK) mktrack.o mkviewer.mod2.o ratint.f Memory2.o -g
	

misc:
	rm -f $(VIEWPOW) 
	rm -f $(OBJVIEW) 
	rm -f $(MKPFILE) 
#	rm -f $(MKTRACK) 
	$(CC) -c getsimparm.c -O1 -mcmodel=medium  -lm $(CAMBOPTION)
	$(CC) -o $(MKPFILE) mkpfile.c header.c -O1 -mcmodel=medium -lm -qopenmp $(CAMBOPTION)
	$(F90C) -c $(MKPOWER).f90 -I$(CAMBDIR) -I$(CAMBDIR2) -g  -mcmodel=medium -qopenmp $(CAMBOPTION)
	$(F90C) -o $(MKPOW) $(MKPOWER).o $(MYSUB).o getsimparm.o header.o -I$(CAMBDIR) -L$(CAMBDIRLIB) -I$(CAMBDIR2) -L$(CAMBDIRLIB2) $(CAMBL) $(CAMBOPTION) -O1 -mcmodel=medium   -qopenmp
ifeq ($(CAMBTYPE), CAMB)
	$(F90C) -c $(VIEWPOWER).f90 -I$(CAMBDIR) -I$(CAMBDIR2) -g  -mcmodel=medium -qopenmp
	$(F90C) -o $(VIEWPOW) $(VIEWPOWER).o -I$(CAMBDIR) -L$(CAMBDIRLIB) -I$(CAMBDIR2) -L$(CAMBDIRLIB2) $(CAMBL)  -g -mcmodel=medium -qopenmp
endif

	$(F90C) -o growthtest.exe growthtest.F -I./$(SUBDIRS) -O1  -mcmodel=medium 
#	$(MAKE) mkview

PCORRTEST:
	rm -rf {PCORRTEST}
	$(CC) -c pcorrtest.c   -O1
	$(F90C) -nofor-main -o $(PCORRTEST) pcorrtest.o $(MYSUB).o getsimparm.o header.o -I$(CAMBDIR) -L$(CAMBDIRLIB) -I$(CAMBDIR2) -L$(CAMBDIRLIB2) $(CAMBL)  $(CAMBOPTION) -O1 -mcmodel=medium   -qopenmp



default: all

#$(MKPOW): $(OBJPOW)



include Rules.make
