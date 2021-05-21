linux: compile execute
mac: compile_mac execute

compile:
	g++  -o main.out test.cpp main2.cpp -std=c++11 -fopenmp -larmadillo

compile_mac:
	c++  -o main.out BoltzPar.cpp main_par.cpp -std=c++11 -larmadillo -lomp -Xpreprocessor -fopenmp

execute:
	./main.out
