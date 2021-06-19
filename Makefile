all:	simu

simu:	main.cpp
	c++ -O2 -ffast-math -o simu main.cpp -std=c++17 -lpthread
