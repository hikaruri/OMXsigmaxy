# For gcc(macports) + intelMKL

FFTROOT = /opt/local
LBSROOT = /opt/intel/mkl
CC = mpicc -fopenmp -O3 -I$(LBSROOT)/include -I$(FFTROOT)/include
LIB= -L$(FFTROOT)/lib -lfftw3 -L$(LBSROOT)/lib -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lpthread -lgfortran -lmkl_scalapack_lp64 -lmkl_blacs_mpich_lp64 -lmpi_mpifh

# For MateriApps LIVE!
# CC = mpicc -O3 -fopenmp -I/usr/include
# LIB = -L/usr/lib -llapack -lblas -lgfortran -lmpi_mpifh

OBJS_sigmaxy = sigmaxy.o read_scfout.o
sigmaxy:	$(OBJS_sigmaxy)
	$(CC) $(OBJS_sigmaxy) $(LIB) -lm -o sigmaxy

sigmaxy.o: sigmaxy.c read_scfout.h 
	$(CC) -c sigmaxy.c

clean: 
	rm -f *.o sigmaxy

