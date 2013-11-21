PROG = wells
CC = gcc 
CFLAGS = -Wall 
OBJ = wells.o
$(PROG): $(OBJ)

wells.o: wells.c wells.h

clean:
	rm -f $(PROG) $(OBJ)

tar:
	tar -cvzf wells.tgz `hg st -c | awk '{print $$2}'` .hg


