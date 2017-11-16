CC = gcc

TARGET = cgp

all: $(TARGET).c
	$(CC) $(TARGET).c -o $(TARGET)

clean:
	rm -f $(TARGET)