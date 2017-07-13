GOOGLETEST_DIR ?= /opt/googletest/1.8.0
HTSLIB_PREFIX ?= /opt/lib/htslib/1.4.1

CXX=g++
CPPFLAGS=-I$(GOOGLETEST_DIR)/include -I$(HTSLIB_PREFIX)/include
CXXFLAGS=-std=c++14 -pthread

ifdef GMOCK_STATIC
	LIBS=$(GOOGLETEST_DIR)/lib/libgmock.a $(HTSLIB_PREFIX)/lib/libhts.a
	LDADDS=-lbz2 -lcurl -lcrypto -lz -llzma
else
	LIBS=$(HTSLIB_PREFIX)/lib/libhts.a
	LDADDS=-lgmock -lbz2 -lcurl -lcrypto -lz -llzma
endif

SOURCES = *.cpp

all: runner


runner: $(SOURCES)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(LDFLAGS) -o $@ $(SOURCES) $(LIBS) $(LDADDS)

check: runner
	./runner

clean:
	@rm -f runner

.PHONY: clean check all
