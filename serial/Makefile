#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = mpic++
FLAGS = -fopenmp 
UTILITIES = main.cpp 
OUT = kdtree.x

SRC = main.cpp utility.cpp  

#$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)

all:	
	$(COMPILER) $(UTILITIES) -g -D double_data -o $(OUT) 

int:
	$(COMPILER) $(UTILITIES) -g -D int_data -o $(OUT)

double:
	$(COMPILER) $(UTILITIES) -g -D double_data -o $(OUT)

optimize:
	$(COMPILER) $(UTILITIES) -g -D double_data -o $(OUT) -O3

debug_int:
	$(COMPILER) $(UTILITIES) -g -D int_data -D int_data_DEBUG -o $(OUT)

debug_double:
	$(COMPILER) $(UTILITIES) -g -D double_data -D double_data_DEBUG -o $(OUT)
 
clean:
	rm *.x
