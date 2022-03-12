#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = mpic++
FLAGS = -fopenmp 
UTILITIES = main_mpi.cpp 
OUT = kdtree_mpi.x

SRC = main.cpp utility.cpp  

#$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

output:	

	$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT) -O3

debug:

	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D DEBUG -o $(OUT) -O3       
 
clean:
	rm *.x
