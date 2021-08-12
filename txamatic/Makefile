# Makefile for txamatic

LIB = -lm -lik
LINC = -Lik
INC = -Iik
CFLAGS = -O2
OBJECTS =

APP = txamatic
OBJ = txamatic.o

default:
	make $(APP)

$(APP): $(OBJ) $(OBJECTS)
	$(CC) -o $(APP) $(OBJ) $(OBJECTS) $(LINC) $(LIB)

clean:
	rm -f *.o $(APP)

%.o: %.c
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
