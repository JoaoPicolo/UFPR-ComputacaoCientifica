    CC     = gcc -g -std=c11 -std=gnu99
    CFLAGS = -Wall
    LFLAGS = -lm

      PROG = matrizInv 
      OBJS = utils.o \
             comunicacaoIO.o \
             sistemasLineares.o \
             matrizes.o

.PHONY: limpa faxina clean purge all

all: matrizInv

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

faxina purge:
	@rm -f *~ *.bak

limpa clean:   faxina
	@rm -f *.o core a.out
	@rm -f $(PROG)

