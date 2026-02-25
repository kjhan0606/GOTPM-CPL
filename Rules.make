#
# Makefile for PM/Tree code namu.exe - Aug 31, 2017 - JHK
#


RANLIB = ranlib



#####################################
# System specific definitions
#####################################

#########################################################
####      KIAS PG Compiler                          #####
#########################################################
AR = ar rcv
FC = mpiifx 
CC = mpiicx 
#F90C = pgf90
F90C = ifx -fpp
##FFTWDIR = /user/kjhan/fftwfinal/
#FFTWDIR = /home/kjhan/fftw-intel/
FFTWDIR = /home/kjhan/local/
#FFTWDIR = /home/mhee7173_snu/local/
##
#MKL_FFTWDIR = /home01/g053pcb/fftw_mkl_intel_openMPI_single/
#MKL_LIB = /applic/compilers/intel/11.1/mkl/lib/em64t

#OPT = -DPGCC -mcmodel=medium -fast -fastsse
#OPT = -DINTEL -mcmodel=medium -O1 -openmp
OPT = -DINTEL -O2 -qopenmp
NVCC = nvcc --maxrregcount 128
#F90C = ifort
#FFTWDIR = /home/t075pcb/intel/
#OPT = -DINTEL -openmp -O2 -xSSE4.2 -m64

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

#####
#####
F90INCLUDES = -I$(CAMBDIR) -I$(CAMBDIR2)
COMFLAGS = -DINDEX -DVarPM   -DXYZDBL  -DMPI_SLOW $(CAMBOPTION)  # -DDEMODEL=$(CAMBDE)

FDFLAGS =  -DINCLUDE_TREE_FORCE  

CUFLAGS = $(CUOPT) -I$(CUDA)/include -L$(CUDA)/lib64/ -I$(CUDA_SDK)/common/inc -L$(CUDA_SDK)/lib

# CUDA OPTIMIZATION
CUOPT = -O3  -use_fast_math  $(COMFLAGS)
# CUDA DIRECTORY
CUDA = /usr/local/cuda/
# CUDA SDK DIRECTORY
CUDA_SDK = /usr/local/cuda/C
CAMBLIBS = -L$(CAMBDIRLIB) -L$(CAMBDIRLIB2) $(CAMBL)  # -pgf90libs


LIBS = -L$(FFTWDIR)/lib/ -lfftw3f_threads -lfftw3f_mpi  -lfftw3f -lfftw3f_omp    $(CAMBLIBS)  # -L$(CUDA)/lib64 -lcudart
# -lsrfftw_mpi -lsfftw_mpi -lsrfftw -lsfftw   $(CAMBLIBS)  # -L$(CUDA)/lib64 -lcudart
#	$(MKL_LIB)/libmkl_intel_lp64.a -Wl,--start-group -L$(MKL_LIB) -lmkl_cdft_core -lmkl_blacs_openmpi_lp64\
#	-lmkl_sequential  -lmkl_core -lmkl_intel_lp64 -Wl,--end-group





#CDFLAGS = -DWGROUPSIZE=500 -DNMEG=21500L -DINCLUDE_TREE_FORCE \

CDFLAGS = -DWGROUPSIZE=80 -DNMEG=8000L -DINCLUDE_TREE_FORCE \
        -D_LARGE_FILES -DSAVESLICE  -DPMSEEDFORCE  # -DUSE_GPU\
		#-DTSC_OLD -DOLD_FDA4


FFLAGS = $(FDFLAGS) $(OPT) $(COMFLAGS)  
CFLAGS = $(OPT) $(CDFLAGS)  $(COMFLAGS) 
LDFLAGS = $(OPT)  -nofor-main
#LDFLAGS = $(OPT)  
#CAMBLIBS = -L./$(CAMBDIR) -lcamb



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
.SUFFIXES: .exe .o .c .f .F .f90 .cu

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
