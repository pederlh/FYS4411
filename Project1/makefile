ser: compile execute
par: compile_par execute

compile:
	c++ -o main.out Solver.cpp Psi.cpp main.cpp -larmadillo -std=c++11

compile_par:
	g++ -o main.out Solver.cpp Psi.cpp main.cpp -std=c++11 -fopenmp

execute:
	./main.out
