CC = gcc
CC_OPTIONS =   -O2 -Wall 
EXEC = euler.x
LIBS = -lm 

OBJS = estructuras.o io.o main.o lax_wendoff.o
INCL = estructuras.h io.h lax_wendoff.h inc.h
OUTS = shock.dat


all: $(OUTS)

.c.o:
	$(CC) $(CC_OPTIONS) -c $<

$(OUTS):$(EXEC)
	./$(EXEC)

$(EXEC):$(OBJS)
	$(CC) $(OBJS) $(LIBS) -o $(EXEC)

clean:
	rm -f $(OBJS) *~ core* ${EXEC}
