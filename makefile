CC = gcc
CC_OPTIONS =   -O2 -Wall 
EXEC = euler.x
LIBS = -lm 

OBJS = main.o estructuras.o io.o lax_wendoff.o struct.o init.o io_sed.o fvm.o space.o 
INTER = sedov.dat tiempo.dat shock.dat
GRAPH = graph.py

all: $(OBJS)

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

plotsedov:
	python $(GRAPH) sedov

plotshock:
	python $(GRAPH) shock

sedov:
	./$(EXEC) sedov

shock:
	./$(EXEC) shock

exec: $(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -f $(OBJS) *~ core* ${EXEC} ${INTER}
