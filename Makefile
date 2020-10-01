.DEFAULT_GOAL := all
./PHONY: clean dir install uninstall count dev_install

INSTALLDIR=/usr/local/bin

CC=gcc
OBJ=./obj
SRC=./src
BIN=./bin
TEST=./test
CFLAGS=-I ./src #-std=c99
LIBS=-lm

all: dir ampd colextract rowextract ampdpreproc

$(OBJ)/%.o: $(SRC)/%.c
	$(CC) -c $(CFLAGS) $< -o $@

ampd: $(OBJ)/ampd.o $(OBJ)/ampdr.o $(OBJ)/filters.o
	$(CC) -o $(BIN)/ampd $(OBJ)/ampd.o $(OBJ)/ampdr.o $(OBJ)/filters.o $(LIBS)

colextract: $(OBJ)/colextract.o
	$(CC) -o $(BIN)/colextract $(OBJ)/colextract.o $(LIBS)

rowextract: $(OBJ)/rowextract.o
	$(CC) -o $(BIN)/rowextract $(OBJ)/rowextract.o $(LIBS)

ampdpreproc: $(OBJ)/ampdpreproc.o $(OBJ)/filters.o
	$(CC) -o $(BIN)/ampdpreproc $(OBJ)/ampdpreproc.o $(OBJ)/filters.o $(LIBS)

dir: 
	mkdir -p $(OBJ)
	mkdir -p $(BIN)

clean:
	rm -f $(OBJ)/*
	rm -f $(BIN)/*
	rm -rf $(TEST)/out/*

count:
	find . -name '*.c' | xargs wc -l
	find . -name '*.h' | xargs wc -l

install:
	cp $(BIN)/ampd $(INSTALLDIR)/ampd
	@echo ''
	@echo 'Installed ampd in /usr/local/bin. Please copy ampdpreproc, colextract, rowextract and scripts into PATH manually if needed.'

uninstall:
	rm $(INSTALLDIR)/ampd

dev_install:
	mkdir -p $${HOME}/bin
	cp $(BIN)/ampd $${HOME}/bin/ampd
	cp $(BIN)/colextract $${HOME}/bin/colextract
	cp $(BIN)/rowextract $${HOME}/bin/rowextract
	@echo 'Make sure ${HOME}/bin is in PATH'

