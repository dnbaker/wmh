CXX?=g++
all: test
%: %.cpp $(wildcard *.h)
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra
