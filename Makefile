CC=gcc
OBJ=./obj
SRC=./src
BIN=./bin
CFLAGS=-I ./src
LIBS=-lm

all: dir ampd colextract

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

ampd: $(OBJ)/ampd.o
	$(CC) -o $(BIN)/ampd $(OBJ)/ampd.o $(LIBS)

colextract: $(OBJ)/colextract.o
	$(CC) -o $(BIN)/colextract $(OBJ)/colextract.o $(LIBS)


dir: 
	mkdir -p $(OBJ)
	mkdir -p $(BIN)

clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*

count:
	find . -name '*.c' | xargs wc -l
	find . -name '*.h' | xargs wc -l
