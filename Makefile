all:	simu

simu:	main.cpp
	c++ -O3 -funroll-loops -fomit-frame-pointer -ffast-math -march=native -o simu  main.cpp -std=c++17 -lpthread
