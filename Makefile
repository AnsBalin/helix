CC=gcc
CFLAGS=
LDFLAGS=-std=c99
SOURCES=src/md_main.c src/md_init.c src/md_algebra.c src/md_forces.c src/md_simulation.c src/md_io.c
OBJECTS=$(SOURCES: .c=.o)
EXECUTABLE=md_out
NAGCDIR=/opt/NAG/clmi623dgl

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