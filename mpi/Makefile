#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = mpic++
FLAGS = -fopenmp
UTILITIES = main_mpi.cpp 
OUT = kdtree.x

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
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D int_data -D DEBUG -o $(OUT)

debug_double:
	$(COMPILER) $(FLAGS) $(UTILITIES) -g -D double_data -D DEBUG -o $(OUT)     
 
clean:
	rm *.x
