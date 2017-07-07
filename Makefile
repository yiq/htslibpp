HTSLIB_PREFIX ?= /opt/lib/htslib/1.4.1

CXXFLAGS=-std=c++14 -I$(HTSLIB_PREFIX)/include
LDFLAGS=-L$(HTSLIB_PREFIX)/lib -Wl,-rpath,$(HTSLIB_PREFIX)/lib
LDADDS=-lhts


all: check

check:
	@make -C test

.PHONY: all check
