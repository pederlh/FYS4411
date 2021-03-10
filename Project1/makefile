ser: compile execute
par: compile_par execute

compile:
	c++ -o main.out Solver.cpp Psi.cpp main_1D.cpp -larmadillo -std=c++11

compile_par:
	g++-10 -o main.out Solver.cpp Psi.cpp main_1D.cpp -std=c++11 -fopenmp

execute:
	./main.out
