CCC = g++       # compiler
AR = ar         # archiver
LD = g++        # linker

ARFLAGS = rcs

LDFLAGS += -Wl,--no-as-needed -fPIC -ggdb3 -O3 -lstdc++
LDFLAGS += $(shell root-config --libs)

#LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaIntern
#LDFLAGS += -L$(COAST_DIR)/lib -lCorsikaFileIO

CXXFLAGS += -c -fPIC  -ggdb3 -O3
#CXXFLAGS += -I$(COAST_DIR)/include
CXXFLAGS += $(shell root-config --cflags)

SRCFILES = $(wildcard *.cc)
OBJECTS = $(patsubst %.cc, %.o, ${SRCFILES})

EXE = fafo

all: ${EXE}

fafo: fafo.o

clean:
	@rm -f *.o *.a *.so *~ core $(EXE)

