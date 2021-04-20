linux: compile execute
mac: compile_mac execute

compile:
	g++  -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

compile_mac:
	c++  -o main.out Solver.cpp Psi.cpp BoltzmannMachine.cpp main.cpp -std=c++11 -larmadillo -lomp -Xpreprocessor -fopenmp

execute:
	./main.out
