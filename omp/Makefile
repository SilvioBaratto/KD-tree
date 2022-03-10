#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = mpiCC
FLAGS = -fopenmp 
UTILITIES = main_omp.cpp 
OUT = kdtree_omp.x

SRC = main.cpp utility.cpp  

#$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

all:	
	$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)          
 
clean:
	rm *.x
