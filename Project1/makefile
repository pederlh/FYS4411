linux: compile execute
mac: compile_mac execute

compile:
	g++  -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

compile_mac:
	g++-10 -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

execute:
	./main.out
