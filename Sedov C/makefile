CC = gcc
CC_OPTIONS =   -O2 -Wall 
EXEC = euler.x
LIBS = -lm 

OBJS = struct.o init.o io.o main.o fvm.o space.o
INCL = Makefile struct.h io.h init.h fvm.h space.h


all: $(EXEC)

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

plotsedov:

sedov:
	./$(EXEC)

$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)
clean:
	rm -f $(OBJS) *~ core* ${EXEC}
