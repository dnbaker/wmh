CXX?=g++
all: bmh
%: %.cpp %.h
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra
