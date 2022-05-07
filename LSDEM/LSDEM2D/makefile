CC = mpicxx 
CC_INCLUDE = -I eigen-3.3.7/Eigen/ 
TARGETS = Main2


all: $(TARGETS)

Main2: maintest.cpp $(wildcard *.h)
	$(CC) $< -Wall -w -O3 -DEIGEN_NO_DEBUG -o $@ -O3 -std=c++11 $(CC_INCLUDE) -fopenmp -lgmp -msse2 -pedantic

clean:
	rm -f $(TARGETS)
	
again: clean $(TARGETS)

