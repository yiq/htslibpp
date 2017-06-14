HTSLIB_PREFIX ?= /opt/lib/htslib/1.4.1
GSL_CONFIG ?= /opt/lib/gsl/2.3/bin/gsl-config

GSL_CFLAGS = $(shell $(GSL_CONFIG) --cflags)
GSL_LIBS   = $(shell $(GSL_CONFIG) --libs-without-cblas)

CXXFLAGS=-std=c++14 -I$(HTSLIB_PREFIX)/include
LDFLAGS=-L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
LDADDS=-lhts


all: main dp_stats random

main: main.cc htslibpp.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDADDS)

random: random.cc htslibpp.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDADDS)

dp_stats: dp_stats.cc htslibpp.h
	$(CXX) $(CXXFLAGS) $(GSL_CFLAGS) $(LDFLAGS) -o $@ $< $(GSL_LIBS) $(LDADDS)


clean:
	rm -f *.o main dp_stats random
