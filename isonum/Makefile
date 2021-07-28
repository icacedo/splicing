# Makefile for isonum

LIB = -lm -lik
LINC = -Lik
INC = -Iik
CFLAGS = -O2 -Wall -Werror
OBJECTS = 

APP = isonum
OBJ = isonum.o

default:
	make $(APP)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LINC) $(LIB)

clean:
	rm -f *.o $(APP)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
