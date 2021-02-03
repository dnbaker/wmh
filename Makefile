CXX?=g++
all: test
%: %.cpp $(wildcard *.h)
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra -march=native

uo%: %.cpp $(wildcard *.h)
	$(CXX) -O1 -std=c++17 -o $@ $< -Wall -Wextra -march=native
p%: %.cpp $(wildcard *.h)
	$(CXX) -O1 -std=c++17 -o $@ $< -Wall -Wextra
