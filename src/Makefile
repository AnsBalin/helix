CC=gcc
CFLAGS=
LDFLAGS=-std=c99
SOURCES=md_test.c md_init.c md_algebra.c md_forces.c md_simulation.c md_io.c
OBJECTS=$(SOURCES: .c=.o)
EXECUTABLE=md_out
NAGCDIR=/home/balin/Library/NAG/cll6i24dcl

all: $(SOURCES) $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $(OBJECTS) -o $@ -I${NAGCDIR}/include \
					-m64 ${NAGCDIR}/lib/libnagc_nag.a \
					${NAGCDIR}/rtl/intel64/libifcore.a ${NAGCDIR}/rtl/intel64/libimf.a \
        			${NAGCDIR}/rtl/intel64/libirc.a \
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