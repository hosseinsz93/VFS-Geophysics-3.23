PETSC_DIR = /home/hossein/lib/petsc-3.23.6/src
PETSC_ARCH = linux-gnu-opt

HYPRE_DIR = /opt/ohpc/pub/libs/gnu7/mvapich2/hypre/2.14.0
OPENBLAS_DIR = /opt/ohpc/pub/libs/gnu7/openblas/0.2.20

TEC360HOME = $(HOME)/lib/TECIO/360_2009_R2
TECINC = $(TEC360HOME)/include
TECLIB = $(TEC360HOME)/lib

CC = mpicxx
EXEC = VFS-Geophysics-3.23

CFLAGS += -I$(PETSC_DIR)/include \
          -I$(PETSC_DIR)/$(PETSC_ARCH)/include \
          -I$(HYPRE_DIR)/include \
          -I$(TECINC)

LDFLAGS += -L$(PETSC_DIR)/$(PETSC_ARCH)/lib \
           -L$(HYPRE_DIR)/lib \
           -L$(OPENBLAS_DIR)/lib

LIBFLAG = -lpetsc -lHYPRE -lopenblas -lpthread -lrt -ldl -lstdc++ -lgfortran

SOURCEC = bcs.c bmv.c canopy.c compgeom.c ibm.c ibm_io.c init.c \
          main.c metrics.c poisson.c rhs.c timeadvancing.c \
          timeadvancing1.c variables.c fsi.c implicitsolver.c \
          fsi_move.c solvers.c rhs2.c wallfunction.c \
          les.c k-omega.c distance.c level.c momentum.c \
          poisson_hypre.c rotor_model.c \
          convection_diffusion.c sediment_transport.c turbinecontrol.c turbinestructure.c wave.c

OBJSC = $(SOURCEC:%.c=%.o)

CPPFLAGS = -DNDEBUG -DTECIO=1 -O3

test: $(OBJSC)
	$(CC) -o $(EXEC) $(OBJSC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

VFS-Geophysics-3.1-mhk: $(OBJSC)
	$(CC) -o VFS-Geophysics-3.23 $(OBJSC) $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

data-3.1: data.o
	$(CC) -o data-3.1 data.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG) $(TECLIB)/tecio64.a

data-small-3.1: data-small.o
	$(CC) -o data-small-3.1 data-small.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG) $(TECLIB)/tecio64.a

data05: data05.o
	$(CC) -o data05 data05.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

itfcsearch: itfcsearch.o
	$(CC) -o itfcsearch itfcsearch.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

hill: hill.o
	$(CC) -o hill hill.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

turbine: TurbineLoc.o
	$(CC) -o TurbineLoc TurbineLoc.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

terrain: terrain.o
	$(CC) -o terrain terrain.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

Tec2UCD: OSLTec2UCD.o
	$(CC) -o Tec2UCD OSLTec2UCD.o $(CPPFLAGS) $(CFLAGS) $(LDFLAGS) $(LIBFLAG)

clean:
	-rm -f *.o
	-rm -f data-3.1
	-rm -f VFS-Geophysics*
