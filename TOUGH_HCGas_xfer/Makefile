# find source code files
FTN_FILES=FTSim_GasEOS_Definitions.f90 FTSim_Allocate_Memory.f90 FTSim_Utility_Functions.f90 FTSim_Gas_Properties.f90 FTSim_Media_Properties.f90 FTSim_Gas_Specifics.f90 FTSim_Main.f90 FTSim_Gas_EOS.f90 FTSim_Matrix_Solvers.f90 FTSim_Executive.f90 FTSim_Inputs.f90
FTN_OBJS=FTSim_GasEOS_Definitions.o FTSim_Allocate_Memory.o FTSim_Utility_Functions.o FTSim_Gas_Properties.o FTSim_Media_Properties.o FTSim_Gas_Specifics.o FTSim_Main.o FTSim_Gas_EOS.o FTSim_Matrix_Solvers.o FTSim_Executive.o FTSim_Inputs.o

#initial flags
FTN_FLAGS=-O3
CXX_FLAGS=-O3 -std=c++11

# specify compiler, assume ICC
ifdef USE_GCC
  FTN=gfortran
  CXX=gcc-6
  FTN_FLAGS+=-cpp -ffree-form -ffree-line-length-0 -Wextra -Wall -pedantic
  CXX_FLAGS+=-Wextra -Wall
else
  FTN=ifort
  CXX=icc
  FTN_FLAGS+=-fpp
endif

# append flags if using OpenMP
ifdef USE_OMP
  ifdef USE_GCC
    FTN_FLAGS+=-fopenmp -DUSE_OMP
    CXX_FLAGS+=-fopenmp -DUSE_OMP
  else
    FTN_FLAGS+=-qopenmp -DUSE_OMP
    CXX_FLAGS+=-qopenmp -DUSE_OMP
  endif
endif

# append flags for C++ and find files
ifdef USE_CPP
  FTN_FLAGS+=-lstdc++ -DUSE_CPP
  CXX_FILES=$(shell ls *.cpp)
  #CXX_OBJS=$(shell ls *.o)
  CXX_OBJS=FTSim_CPP_Funcs.o
endif

# append flags to use timers
ifdef USE_TIMER
  FTN_FLAGS+=-DUSE_TIMER
  CXX_FLAGS+=-DUSE_TIMER
endif

all:
	for file in $(CXX_FILES); do \
		$(CXX) $(CXX_FLAGS) -c $$file; \
	done
	for file in $(FTN_FILES); do \
		$(FTN) $(FTN_FLAGS) -c $$file; \
	done
	$(FTN) $(FTN_FLAGS) -o FTSim.exe $(CXX_OBJS) $(FTN_OBJS)

clean:
	git checkout FracGrid/MESH
	rm *.o FTSim.exe *.mod ALLOC LINEQ SAVE VERS INCON Plot* GENER
