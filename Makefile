# the compiler to use
CC = gcc
# compiler flags: -g adds debugging info, -Wall turns on most warnings
CFLAGS = -g -Wall
# the build target executable name
TARGET = raycell
LIBS = -lraylib -lm

all: $(TARGET)

$(TARGET): raycell.c ca_core.c
	$(CC) $(CFLAGS) -o $(TARGET) raycell.c ca_core.c $(LIBS)

clean:
	rm -f $(TARGET) *.o
