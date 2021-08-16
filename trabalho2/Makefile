# ###################################
# Gabriel Marczuk Thá  - GRR20186070 
# João Pedro Picolo    - GRR20182659
# ###################################

    CC     = gcc -g -std=c11 -std=gnu99
    CFLAGS = -Wall -DLIKWID_PERFMON -I/home/soft/likwid/include -O3 -mavx2 -march=native
    LFLAGS = -lm -L/home/soft/likwid/lib -llikwid

      PROG = calcPolinomial
      OBJS = utils.o \
             tabela.o \
             sistemaLinear.o \
             comunicacaoIO.o

.PHONY: limpa faxina clean purge all

all: calcPolinomial

%.o: %.c %.h
	$(CC) -c $(CFLAGS) $<

$(PROG) : % :  $(OBJS) %.o
	$(CC) -o $@ $^ $(LFLAGS)

faxina purge:
	@rm -f *~ *.bak

limpa clean:   faxina
	@rm -f *.o core a.out
	@rm -f $(PROG)

