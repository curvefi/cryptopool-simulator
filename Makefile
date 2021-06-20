all:	simu

simu:	main.cpp
	c++ -O2 -ffast-math -o simu -march=native main.cpp -std=c++17 -lpthread
