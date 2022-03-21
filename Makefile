#CXX = mpic++
.DEFAULT_GOAL := double
COMPILER = 
FLAGS =  
UTILITIES = 
HPC_FLAGS = -O3 -march=native
OUT = 

# OpenMP or MPI?
ifeq ($(src), mpi)
	SRC = mpi/main_mpi.cpp
	FLAGS += -std=c++17 -fopenmp
	COMPILER = mpic++
	OUT = bin/tree_mpi.x
endif
ifeq ($(src), serial)
	SRC = serial/main.cpp
	COMPILER = c++
	FLAGS = -std=c++17 -fopenmp
	OUT = bin/tree_serial.x
endif
ifeq ($(src), omp)
	SRC = omp/main_omp.cpp
	FLAGS += -std=c++17 -fopenmp
	COMPILER = g++
	OUT = bin/tree_omp.x
endif

int:
	$(COMPILER) $(FLAGS) $(HPC_FLAGS) $(UTILITIES) $(SRC) -g -D int_data -o $(OUT)

double:
	$(COMPILER) $(FLAGS) $(HPC_FLAGS) $(UTILITIES) $(SRC) -g -D double_data -o $(OUT)

debug_int:
	$(COMPILER) $(FLAGS) $(HPC_FLAGS) $(UTILITIES) $(SRC) -g -D int_data -D int_data_DEBUG -o $(OUT)

debug_double:
	$(COMPILER) $(FLAGS) $(HPC_FLAGS) $(UTILITIES) $(SRC) -g -D double_data -D double_data_DEBUG -o $(OUT)        
 
clean:
	rm *.x