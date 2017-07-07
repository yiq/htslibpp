GOOGLETEST_DIR ?= /opt/googletest/1.8.0
HTSLIB_PREFIX ?= /opt/lib/htslib/1.4.1

CPPFLAGS=-I$(GOOGLETEST_DIR)/include -I$(HTSLIB_PREFIX)/include
CXXFLAGS=-std=c++14
LDFLAGS=-L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
LDADDS=-lhts
LIBS=$(GOOGLETEST_DIR)/lib/libgmock.a

SOURCES = runner.cpp \
	  htsfile.cpp \
	  bam_header.cpp

all: check


runner: $(SOURCES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SOURCES) $(LIBS) $(LDADDS)

check: runner
	DYLD_LIBRARY_PATH=/opt/lib/htslib/1.4.1/lib ./runner

clean:
	rm runner

.PHONY: clean check all
