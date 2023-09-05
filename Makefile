cc := gcc
LIB = src/alloc.c src/log.c src/file.c
CFLAGS = -lz -fopenmp -O3 -Dparallel -g
TARGET = roa

all: $(TARGET)

$(TARGET): src/main.c $(LIB)
	$(cc) -o $(TARGET) src/main.c $(LIB) $(CFLAGS)

clean:
	rm -rf $(TARGET)

.PHONY:
	clean
