CXX             = g++
CXXFLAGS_PAR    = -std=c++17 -O3 -fopenmp
PAR_SRCS        = ParallelBipartiteButterflies.cpp
PAR_TARGET      = ParallelBipartiteButterflies

CXXFLAGS_SEQ    = -std=c++17
SEQ_SRCS        = SequentialBipartiteButterflies.cpp
SEQ_TARGET      = SequentialBipartiteButterflies

all: $(PAR_TARGET) $(SEQ_TARGET)

$(PAR_TARGET): $(PAR_SRCS)
	$(CXX) $(CXXFLAGS_PAR) -o $@ $<

$(SEQ_TARGET): $(SEQ_SRCS)
	$(CXX) $(CXXFLAGS_SEQ) -o $@ $<

clean:
	rm -f $(PAR_TARGET) $(SEQ_TARGET)
