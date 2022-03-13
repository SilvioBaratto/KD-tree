#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = g++
FLAGS = -fopenmp 
UTILITIES = main_omp.cpp 
OUT = kdtree_omp.x

SRC = main.cpp utility.cpp  

#$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

all:	
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D double_data -o $(OUT) 

int:
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D int_data -o $(OUT)

double:
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D double_data -o $(OUT)

optimize:
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D double_data -o $(OUT) -O3

debug_int:
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D int_data -D int_data_DEBUG -o $(OUT)

debug_double:
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D double_data -D double_data_DEBUG -o $(OUT)        
 
clean:
	rm *.x
