# Compiler and flags
CC = gcc
CFLAGS = -O3 -fopenmp -Wall
TARGET = spmv

# Source files
SRC = draft1.c specifications.c mmio.c
OBJ = $(SRC:.c=.o)

# Default target
all: $(TARGET)

# Link object files into executable
$(TARGET): $(OBJ)
	$(CC) $(CFLAGS) -o $@ $(OBJ)

# Compile each .c file into .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up
clean:
	rm -f $(OBJ) $(TARGET)

# Run the program with a .mtx file
run: $(TARGET)
	OMP_NUM_THREADS=4 ./$(TARGET) matrix.mtx