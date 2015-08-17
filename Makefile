CC=gcc
# push test
CFLAGS=
LDFLAGS=-std=c99
SOURCES=md_test.c md_init.c md_algebra.c md_forces.c md_simulation.c md_io.c
OBJECTS=$(SOURCES: .c=.o)
EXECUTABLE=md_out

# Comment out or add directory to NAG file:
NAGCDIR=/opt/NAG/clmi623dgl # Andrew Macbook
# NAGCDIR=/home/balin/Library/NAG/cll6i24dcl # Andrew's Linux

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ -I${NAGCDIR}/include \
					-m64 ${NAGCDIR}/lib/libnagc_nag.a \
					-ldl -lpthread -lm

.cpp.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm md_out


#objects = md_test.o md_init.o md_forces.o md_algebra.o 
#
#all: $(objects)
#	$(CC) $(LDFLAGS) $(objects) -o $(EXECUTABLE) -lm
#
#md_test.o: md_test.c md_algebra.c md_forces.c md_init.c
#md_forces.o: md_forces.c md_algebra.c
#md_init.o: md_init.c
#md_algebra.o: md_algebra.c
#
#clean :
#	rm edit $(objects)
