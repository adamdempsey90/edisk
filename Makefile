EXECUTABLE=edisk
SOURCES=init.c alloc.c main.c readinputs.c utils.c output.c boundary.c restart.c history.c algo_driver.c implicit.c transport.c rktvd.c poisson.c
HEADER=edisk.h rk45.h

LDFLAGS=-llapack -lblas -lm -lgomp

CFLAGS=-c -fopenmp -Wall -O3 -g 


BIN=bin/
SRC=src/
IN=inputs/
PY=src/pyutils/


UNAME=`uname`

ifeq ($(UNAME),Linux)
CC=gcc
endif

ifeq ($(UNAME),Darwin)
CC=gcc-4.9
endif

#!!!!!DO NOT EDIT ANYTHING UNDER THIS LINE!!!!!
OBJECTS=$(SOURCES:.c=.o)
CSOURCES=$(addprefix $(SRC),$(SOURCES))
COBJECTS=$(addprefix $(BIN),$(OBJECTS))
CHEADER=$(addprefix $(SRC),$(HEADER))




all: $(CSOURCES) $(EXECUTABLE) 

$(EXECUTABLE): $(COBJECTS) 
	$(CC)  $(COBJECTS) $(LDFLAGS) -o $@

$(BIN)%.o: $(SRC)%.c $(CHEADER) 
	$(CC) $(CFLAGS) $< -o $@

	
clean:
	rm $(COBJECTS) $(EXECUTABLE) 
