HTSLIB_PREFIX ?= /opt/htslib/1.3.1

CXXFLAGS=-std=c++14 -I$(HTSLIB_PREFIX)/include
LDFLAGS=-L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
LDADDS=-lhts


all: main dp_stats

main: main.cc htslibpp.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDADDS)

dp_stats: dp_stats.cc htslibpp.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDADDS)


clean:
	rm -f *.o main
