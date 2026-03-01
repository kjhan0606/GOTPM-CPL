#
# Makefile for PM/Tree code pyul.exe - April 23, 2001 - JJD
#


RANLIB = ranlib



#####################################
# System specific definitions
#####################################

#########################################################
####      KIAS PG Compiler                          #####
#########################################################
AR = ar rcv
FC = mpif90
CC = mpicc
F90C = pgf90
#FFTWDIR = /user/kjhan/fftwfinal/
FFTWDIR = /user/kjhan/fftw-intel/
OPT = -g -DINTEL # -mcmodel=medium -Mlarge_arrays -tp k8-64 -O3 -Mnontemporal # -mp

#########################################################
####      KISTI IBM SP MACHINES                     #####
#########################################################
#AR = ar -X64 rcv
#FC = mpxlf_r
#CC = mpxlc_r
#F90C = mpxlf_r
#FFTWDIR = /user/kjhan/fftwfinal/
#OPT = -O3  -q64 -qtune=pwr5 -qarch=pwr5

#########################################################
####      KISTI IBM SP MACHINES                     #####
#########################################################
#AR = ar rcv
#FC = mpif77
#CC = mpicc
#F90C = ifort
#FFTWDIR = /user/kjhan/fftwfinal/
#OPT = -O1 -xP -ip


#############################################
# List of compilation directives
#############################################
# 
# -DNMEG=nnn      - number of megabytes of storage per processor - default=40
# -DWGROUPSIZE=nnn      - size of the subgroup for sequential I/O of data
# -DINCLUDE_TREE_FORCE - define if you want tree code corrections
# -DTREE               - define if you want tree code corrections
# -DDEBUG              - output debugging information
# -DBIT64              - use on machines using 64 bit addressing 
#                       - i.e. Compaq alpha
# -DAIX                - use on IBM SP3 machines
# -DOBSERVATION_MODE   - save halo data at every ObservationStride'th step
# -DSAVESLICE          - save slice density file at every time step
#############################################################

INCLUDES = -I$(FFTWDIR)/include

#LIBS = -L$(FFTWDIR)/lib/ -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw   $(CAMBLIBS) 
#####
#####
F90INCLUDES = -I./$(CAMBDIR)/
COMFLAGS = -DINDEX -DVarPM  -DXYZDBL 

FDFLAGS =  -DINCLUDE_TREE_FORCE  

CDFLAGS = -DWGROUPSIZE=100 -DNMEG=3000L -DINCLUDE_TREE_FORCE \
        -D_LARGE_FILES -DSAVESLICE  -DPMSEEDFORCE \
		#-DTSC_OLD -DOLD_FDA4


FFLAGS = $(FDFLAGS) $(OPT) $(COMFLAGS)  
CFLAGS = $(OPT) $(CDFLAGS)  $(COMFLAGS) 
LDFLAGS = $(OPT)  -nofor-main
CAMBLIBS =



#################################
# Compaq Alpha
#################################
#FC = f77
#CC = cc 
#INCLUDES = -I/home/kjhan/dolphin/fftw/include
#FDFLAGS = -DBIT64 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT
#DFLAGS = -DBIT64 -DNMEG=256 -DINCLUDE_TREE_FORCE -DTREE -DCOMPACT -DTREEFIX
#FFLAGS = $(FDFLAGS)  -fast -nofor_main
#CFLAGS = $(DFLAGS)  -fast
#LDFLAGS =  -fast -nofor_main
#LIBS = -L/home/kjhan/dolphin/fftw/lib -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw -lm -lz -lmpi -lm
##################################

#--- C Compiler information
#  Leave the rest untouched


#--- Suffix-based compilation rules
.SUFFIXES: .exe .o .c .f .F .f90

#rules to build binary from source


.c.o :
	$(CC) $(CFLAGS) $(INCLUDES) -c $<

.f90.o :
	$(F90C) $(FFLAGS) $(INCLUDES) $(F90INCLUDES) -c $<

.f.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.for.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.F.o :
	$(FC) $(FFLAGS) $(INCLUDES) -c $<

.cu.o :
	$(NVCC) $(CUFLAGS) $(INCLUDES) -c $<

.c.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

.f.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.f90.exe :
	$(F90C) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.for.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.F.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS) $(INCLUDES)

.cu.exe :
	$(FC) $(LDFLAGS) $(INCLUDES) -o $*.exe $(OBJS) $(MPI_LINKS) $(LIBS)

#--- Targets
