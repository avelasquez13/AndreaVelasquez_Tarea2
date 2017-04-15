CC = gcc
CC_OPTIONS =   -O2 -Wall 
EXEC_sed = sedov.x
LIBS = -lm 

OBJS_sed = struct.o init.o io.o main.o fvm.o space.o
INCL_sed = struct.h io.h init.h fvm.h space.h
INTER_sed = sedov.dat tiempo.dat
GRAPH_sed = graphsedov.py
OUTS_sed = prom_rho.png

all: $(OBJS_sed)

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

plotsedov:
	python $(GRAPH_sed)

sedov:
	./$(EXEC_sed)

exec: $(OBJS_sed)
	$(CC) $(OBJS_sed) $(LIBS) -o $(EXEC_sed)

clean:
	rm -f $(OBJS_sed) *~ core* ${EXEC_sed} data/${INTER_sed}
