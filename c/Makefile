CC       := gcc
CFLAGS   := -Wall -g -fsanitize=address -DDEBUG
OBJ      := matrix.o

all: $(OBJ)

release: CFLAGS := -Wall -O2 -DNDEBUG
release: $(OBJ)

%.o: %.c
	$(CC) $(CFLAGS) -c -o $@ $<

clean:
	rm *.o
