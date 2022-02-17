#CXX = mpic++
.DEFAULT_GOAL := output
COMPILER = c++
FLAGS = -fopenmp 
UTILITIES = main.cpp utility.cpp
OUT = kdtree.x

SRC = main.cpp utility.cpp  

all:
	$(COMPILER) $(FLAGS) $(UTILITIES) -o $(OUT)          
 
clean:
	rm *.x
