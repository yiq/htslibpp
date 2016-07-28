CXX=clang++-3.5
#CXX=g++-5
CXXFLAGS=-std=c++14 -stdlib=libc++ -Ihtslib-1.3.1
#CXXFLAGS=-std=c++14 -Ihtslib-1.3.1
LDFLAGS=-Lhtslib-1.3.1 -Wl,-rpath,$(shell pwd)/htslib-1.3.1
LDADDS=-lhts

all: main

main: main.cc htslibpp.h
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $@ $< $(LDADDS)

clean:
	rm -f *.o main
