all:	simu

simu:	main.cpp
	c++ -O3 -funroll-loops -fomit-frame-pointer -ffast-math -march=native -o simu  main.cpp -std=c++17 -lpthread

simud:	main.cpp
	c++ -g -o simud  -fsanitize=address main.cpp -std=c++17 -lpthread

simus:	main-stableswap.cpp
	c++ -O3 -funroll-loops -fomit-frame-pointer -ffast-math -march=native -o simus  main-stableswap.cpp -std=c++17 -lpthread

simudyn:	main-stableswap-dynrate.cpp
	c++ -O3 -funroll-loops -fomit-frame-pointer -ffast-math -march=native -o simudyn  main-stableswap-dynrate.cpp -std=c++17 -lpthread

simul:	main-legacy.cpp
	c++ -O3 -funroll-loops -fomit-frame-pointer -ffast-math -march=native -o simul  main-legacy.cpp -std=c++17 -lpthread
