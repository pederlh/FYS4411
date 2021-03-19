ser: compile execute
par: compile_par execute
par_mac: compile_par_mac execute

compile:
	c++ -o main.out Solver.cpp Psi.cpp main.cpp  -std=c++11

compile_par:
	g++  -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

compile_par_mac:
	g++-10 -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

execute:
	./main.out
