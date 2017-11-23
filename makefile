CC = gcc

FLAGS = -lm

all:
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)

run:
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)
	./main

clean:
	rm -f main