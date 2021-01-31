CXX?=g++
all: test
%: %.cpp $(wildcard *.h)
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra -march=native -Wno-ignored-qualifiers
rel%: %.cpp $(wildcard *.h)
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra -march=native -Wno-ignored-qualifiers -DNDEBUG
m%: %.cpp $(wildcard *.h)
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra -march=native -Wno-ignored-qualifiers -fsanitize=address

uo%: %.cpp $(wildcard *.h)
	$(CXX) -O1 -std=c++17 -o $@ $< -Wall -Wextra -march=native
p%: %.cpp $(wildcard *.h)
	$(CXX) -O1 -std=c++17 -o $@ $< -Wall -Wextra
m%: %.cpp $(wildcard *.h)
	$(CXX) -O3 -std=c++17 -o $@ $< -Wall -Wextra -march=native -Wno-ignored-qualifiers -fsanitize=address

uo%: %.cpp $(wildcard *.h)
	$(CXX) -O1 -std=c++17 -o $@ $< -Wall -Wextra -march=native
p%: %.cpp $(wildcard *.h)
	$(CXX) -O1 -std=c++17 -o $@ $< -Wall -Wextra
