#makefile

DEBUG ?= 0
STATIC ?= 0

# Submodules
PWD = $(shell pwd)
EBROOTHTSLIB ?= ${PWD}/src/htslib/
SDSL_ROOT ?= ${PWD}/src/sdsl-lite/src/sdsl-lite

# Install dir
prefix = ${PWD}
exec_prefix = $(prefix)
bindir ?= $(exec_prefix)/bin


CXX=g++
CXXFLAGS += -std=c++11 -isystem ${EBROOTHTSLIB} -I ${SDSL_ROOT}/include -pedantic -W -Wall -fvisibility=hidden
LDFLAGS += -L ${SDSL_ROOT}/lib -lsdsl -ldivsufsort -ldivsufsort64 -lboost_iostreams -lboost_filesystem -lboost_system -lboost_program_options -lboost_date_time -L${EBROOTHTSLIB} -L${EBROOTHTSLIB}/lib -lpthread



# Additional flags for release/debug
ifeq (${STATIC}, 1)
	LDFLAGS += -static -static-libgcc -pthread -lhts -lz -llzma -lbz2
else
	LDFLAGS += -lhts -lz -llzma -lbz2 -Wl,-rpath,${EBROOTHTSLIB}
endif
ifeq (${DEBUG}, 1)
	CXXFLAGS += -g -O0 -fno-inline -DDEBUG
else ifeq (${DEBUG}, 2)
	CXXFLAGS += -g -O0 -fno-inline -DPROFILE
	LDFLAGS += -lprofiler -ltcmalloc
else
	CXXFLAGS += -O3 -fno-tree-vectorize -DNDEBUG
endif
ifeq (${EBROOTHTSLIB}, ${PWD}/src/htslib/)
	SUBMODULES += .htslib
endif
ifeq (${SDSL_ROOT}, ${PWD}/src/sdsl-lite/)
	SUBMODULES += .sdsl
endif


# External sources
SDSLSOURCES = $(wildcard src/sdsl-lite/lib/*.cpp)
HTSLIBSOURCES = $(wildcard src/htslib/*.c) $(wildcard src/htslib/*.h)
SOURCES = $(wildcard src/*.h) $(wildcard src/*.cpp)

# Targets
BUILT_PROGRAMS = src/validator
TARGETS = ${SUBMODULES} ${BUILT_PROGRAMS}

all:   	$(TARGETS)

.sdsl: $(SDSLSOURCES)
	if [ -r src/sdsl-lite/install.sh ]; then cd src/sdsl-lite/ && ./install.sh ${PBASE}src/sdsl-lite && cd ../../ && touch .sdsl; fi
	if [ ! -r src/sdsl-lite/install.sh ]; then cd src/ && git clone https://github.com/simongog/sdsl-lite.git && cd sdsl-lite/ && ./install.sh ${PBASE}src/sdsl-lite && cd ../../ && touch .sdsl; fi

.htslib: $(HTSLIBSOURCES)
	if [ -r src/htslib/Makefile ]; then cd src/htslib && make && make lib-static && cd ../../ && touch .htslib; fi
	if [ ! -r src/htslib/Makefile ];then cd src/ && git clone https://github.com/samtools/htslib.git && cd htslib && make && make lib-static && cd ../../ && touch .htslib; fi

src/validator: ${SUBMODULES} $(SOURCES)
	$(CXX) $(CXXFLAGS) $@.cpp -o $@ $(LDFLAGS)

install: ${BUILT_PROGRAMS}
	mkdir -p ${bindir}
	install -p ${BUILT_PROGRAMS} ${bindir}

clean:
	if [ -r src/htslib/Makefile ]; then cd src/htslib && make clean; fi
	if [ -r src/sdsl-lite/install.sh ]; then cd src/sdsl-lite/ && ./uninstall.sh && cd ../../ && rm -rf src/sdslLite/; fi
	rm -f $(TARGETS) $(TARGETS:=.o) ${SUBMODULES}

distclean: clean
	rm -f ${BUILT_PROGRAMS}

.PHONY: clean distclean install all
