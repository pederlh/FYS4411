linux: compile execute
mac: compile_mac execute

compile:
	g++  -o main.out BoltzmannMachine.cpp main.cpp -std=c++11 -fopenmp -larmadillo

compile_mac:
	c++  -o main.out BoltzmannMachine.cpp main.cpp -std=c++11 -larmadillo -lomp -Xpreprocessor -fopenmp

execute:
	./main.out
