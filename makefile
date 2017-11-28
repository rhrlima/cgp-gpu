CC = gcc

FLAGS = -lm

all:
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)

run:
	$(CC) src/cgp.c src/main.c -o main $(FLAGS)
	./main

test:
	$(CC) src/cgp.c src/test1.c -o test1
	$(CC) src/cgp.c src/test2.c -o test2
	./test1
	./test2

clean:
	rm -f main