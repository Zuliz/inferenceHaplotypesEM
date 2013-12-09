SRC = src/
INC = inc/
EXE = exe/

CC = gcc
DEBUG = 1
ifeq ($(DEBUG), 1)
    FLAGS = -Wall -pedantic -ansi -g
else
    FLAGS = -Wall -pedantic -ansi
endif
EXEC = $(EXE)main $(EXE)inference_main $(EXE)verification


all: $(EXE) $(EXEC)

$(EXE)verification: $(SRC)verification.c $(SRC)fonctions_verification.c
	@$(CC) $^ -o $@ $(FLAGS)

$(EXE)inference_main: $(SRC)inference_main.c $(SRC)generation_haplotypes.c $(SRC)liste_doublement_chainee.c $(SRC)inference_haplotype_EM.c
	@$(CC) $^ -lm -o $@ $(FLAGS)

$(EXE)main: $(SRC)fonctions.c $(SRC)main.c
	@$(CC) $^ -o $@ $(FLAGS)

$(EXE):
	@mkdir $@

.PHONY: clean

clean:
	@rm -rf $(EXEC)
