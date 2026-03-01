#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#

#--- Executable

PROGS =  opfof.exe 

#--- Objects to build

####################################################################################
# first version
#OBJS		= opfof.omp.o Time.o Treewalk.fof.ordered.o  read.o  header.o
#CDFLAGS = -DWGROUPSIZE=10 -DNMEG=48000L 
####################################################################################
# seoncd version
OBJS		= opfof.omp4.o Time.o Treewalk.fof.ordered.master.o  read.omp.o  header.o fof.eunha.o fof.old.o 
# opfof.omp4.c: experimental & not working now


###################################################################################
#CDFLAGS = -DMulOut -DmaxReads=5 -DOLD_OMP 
CDFLAGS = -DmaxReads=5  -DOLD_OMP #-DMASTERWRITE
####################################################################################
##      COMPILATION OPTIONS
####################################################################################
# -DMulOut : each slave saves the halo data to its own file.
# -DmaxReads : the maximum number of synchronous reading
# -DMASTERWRITE : the master node takes the job saving the halo data (data transfer from slaves)
# -DOLD_OMP : Original Tree Builder. If not selected, the Eunha code OMP-version of 
#          tree Builder (fof.eunha.c; EXPERIMENTAL & NOT WORKING NOW) is used. 
# If the serial writing of the system is slow, choose -DMulOut.
# else if the file-locking in writing is not guaranteed (concurrency issue), choose -DMASTERWRITE.
####################################################################################

####################################################################################
##      PARAMETERS
####################################################################################
# MaxLinkedParticles in "fof.h":  important to set the maximum number of FoF members.


OS = LINUX_GCC



#--- Targets

CAMBDIR = 
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


#mkpower:
#	rm -f $(MKPOW) 
#	rm -f $(OBJPOW) 
#	$(MAKE) $(MKPOW)

default: all

#$(MKPOW): $(OBJPOW)



include Rules.make
