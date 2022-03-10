#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = mpiCC
FLAGS = -fopenmp 
UTILITIES = main_mpi.cpp 
OUT = kdtree_mpi.x

SRC = main.cpp utility.cpp  

#$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

output:	

	$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

debug:

	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D DEBUG -o $(OUT)          
 
clean:
	rm *.x
