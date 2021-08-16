CC = gcc
CFLAGS = -g -Wall
OBJS = operations.o manipulation.o analise.o
DEPS = operations.h manipulation.h types.h

all: analise

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

# Cria o executavel analise
analise: $(OBJS)
	$(CC) -o $@ $^ $(CFLAGS)

clean:
	@rm -f *.o

purge: clean
	@rm analise

# EOF ----------------------------------------------------------%