# Compiler and flags
CC := g++
CFLAGS := -Wall -Wextra -I. -std=c++11 -march=native

# Files
SRCS := main.cpp whfastfpga.cpp whfastfpga_kernel.cpp whfast512.cpp whfast512_kernel.cpp util.cpp
OBJS := $(SRCS:.c=.o)
TARGET := main

# Default target
all: $(TARGET)

# Link objects into the final binary
$(TARGET): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# Compile .c to .o
%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

# Clean up build artifacts
.PHONY: clean
clean:
	rm -f $(OBJS) $(TARGET)
