CC = gcc
CFLAGS = -std=c99 -Wall

main: main.o
	$(CC) -o $@ $^ -lm

main.o: main.c
	$(CC) -c $^

clean:
	rm *.o main
