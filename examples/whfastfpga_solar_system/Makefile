# Compiler and flags
CC        := g++
CFLAGS    := -Wall -Wextra -std=c++11 -march=native -g
SRCDIR    := src
OBJDIR    := build

# Tell the compiler to look in src/ for headers
INCLUDES  := -I$(SRCDIR)

# Source files (just the basenames)
LIB_SRCS  := whfastfpga.cpp \
             whfastfpga_kernel.cpp \
             whfast512.cpp \
             whfast512_kernel.cpp \
             util.cpp

MAIN_SRC  := main.cpp
TEST_SRC  := test.cpp

# Object files go into build/, mirror the basename of each .cpp
MAIN_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAIN_SRC)  $(LIB_SRCS))
TEST_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(TEST_SRC)  $(LIB_SRCS))

# Final executables
TARGET       := main
TEST_TARGET  := test

# Default target: build 'main'
.PHONY: all
all: $(TARGET)

# Link 'main'
$(TARGET): $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

# 'make test' builds the 'test' executable
.PHONY: test
test: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

# Ensure build/ exists, then compile any src/foo.cpp → build/foo.o
$(OBJDIR):
	mkdir -p $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean out all build artifacts
.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(TARGET) $(TEST_TARGET)