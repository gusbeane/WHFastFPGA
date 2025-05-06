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
			 whfast512_constants.cpp \
             util.cpp

MAIN_SRC  := main.cpp
TEST_SRC  := test.cpp
GENGOLD_SRC := gen_golden.cpp
TESTGOLD_SRC := test_golden.cpp

# Object files go into build/, mirror the basename of each .cpp
MAIN_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAIN_SRC)  $(LIB_SRCS))
TEST_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(TEST_SRC)  $(LIB_SRCS))
GENGOLD_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(GENGOLD_SRC)  $(LIB_SRCS))
TESTGOLD_OBJS := $(patsubst %.cpp,$(OBJDIR)/%.o,$(TESTGOLD_SRC)  $(LIB_SRCS))

# Final executables
TARGET       := main
TEST_TARGET  := test
GENGOLD_TARGET := gen_golden
TESTGOLD_TARGET := test_golden

# Default target: build 'main'
.PHONY: all
all: $(TARGET)

# Link 'main'
$(TARGET): $(MAIN_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

# 'make test' builds the 'test' executable
# Overwrite 'make test' to build and run 'test_golden'
.PHONY: test
test: $(TESTGOLD_TARGET)

$(TESTGOLD_TARGET): $(TESTGOLD_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

# Add gen_golden target
.PHONY: gengold
gengold: $(GENGOLD_TARGET)

$(GENGOLD_TARGET): $(GENGOLD_OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $^

# Ensure build/ exists, then compile any src/foo.cpp → build/foo.o
$(OBJDIR):
	mkdir -p $@

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

# Clean out all build artifacts
.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(TARGET) $(TEST_TARGET) $(GENGOLD_TARGET)