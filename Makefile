CC=gcc
OBJ=./obj
SRC=./src
BIN=./bin
CFLAGS=-I ./src
LIBS=-lm

all: dir ampd colextract ampdpreproc

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

ampd: $(OBJ)/ampd.o $(OBJ)/ampdr.o $(OBJ)/filters.o
	$(CC) -o $(BIN)/ampd $(OBJ)/ampd.o $(OBJ)/ampdr.o $(OBJ)/filters.o $(LIBS)

colextract: $(OBJ)/colextract.o
	$(CC) -o $(BIN)/colextract $(OBJ)/colextract.o $(LIBS)

ampdpreproc: $(OBJ)/ampdpreproc.o
	$(CC) -o $(BIN)/ampdpreproc $(OBJ)/ampdpreproc.o $(OBJ)/filter.o $(LIBS)

dir: 
	mkdir -p $(OBJ)
	mkdir -p $(BIN)

clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*

count:
	find . -name '*.c' | xargs wc -l
	find . -name '*.h' | xargs wc -l
