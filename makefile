all:
	g++ -o stochastic -Wall -O2 stochastic\ source/main.cpp stochastic\ source/file-io.cpp
	g++ -o deterministic -Wall -O3 deterministic\ source/main.cpp deterministic\ source/functions.cpp
	g++ -o analysis/ofeatures -Wall -O2 analysis/sources/ofeatures.cpp
	g++ -o analysis/smoothing -Wall -O2 analysis/sources/smoothing.cpp
	g++ -o analysis/t-test -Wall -O2 analysis/sources/t-test.cpp

