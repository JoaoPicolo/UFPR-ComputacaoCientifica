    CC     = gcc -g -std=c11 -std=gnu99
    CFLAGS = -Wall
    LFLAGS = -lm

      PROG = labEqDif 
      OBJS = utils.o \
             EquacoesOrdinarias.o \
             EquacoesParciais.o

.PHONY: limpa faxina clean purge all

%.o: %.c %.h utils.h
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

limpa clean:
	@rm -f *~ *.bak

faxina purge:   limpa
	@rm -f *.o core a.out
	@rm -f $(PROG)

