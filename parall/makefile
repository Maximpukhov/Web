CC=mpicxx
CFLAGS=-O3 -msse2 -funroll-loops
LDFLAGS=

SOURCES=mc_growth.cpp cycle.cpp deform.cpp deform2.cpp deform3.cpp geometry.cpp jump.cpp output.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=mc_simulation

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@

.cpp.o:
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f *.o $(EXECUTABLE)