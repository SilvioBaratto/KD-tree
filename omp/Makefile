#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = mpic++
FLAGS = -fopenmp 
UTILITIES = main_omp.cpp 
OUT = kdtree_omp.x

SRC = main.cpp utility.cpp  

#$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

all:	
	$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT) -O3         
 
clean:
	rm *.x
