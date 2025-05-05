# Compiler and flags
CC      := g++
CFLAGS  := -Wall -Wextra -I. -std=c++11 -march=native -g

# Directories
OBJDIR  := build

# Sources
LIB_SRCS    := whfastfpga.cpp \
               whfastfpga_kernel.cpp \
               whfast512.cpp \
               whfast512_kernel.cpp \
               util.cpp

MAIN_SRC    := main.cpp
TEST_SRC    := test.cpp

# Object lists (in build/)
MAIN_OBJS   := $(patsubst %.cpp,$(OBJDIR)/%.o,$(MAIN_SRC)  $(LIB_SRCS))
TEST_OBJS   := $(patsubst %.cpp,$(OBJDIR)/%.o,$(TEST_SRC)  $(LIB_SRCS))

# Targets
TARGET      := main
TEST_TARGET := test

# Default
.PHONY: all
all: $(TARGET)

# Link main
$(TARGET): $(MAIN_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# Link test
.PHONY: test
test: $(TEST_TARGET)

$(TEST_TARGET): $(TEST_OBJS)
	$(CC) $(CFLAGS) -o $@ $^

# Rule: build directory
$(OBJDIR):
	mkdir -p $@

# Compile any .cpp → build/%.o and ensure build/ exists
$(OBJDIR)/%.o: %.cpp | $(OBJDIR)
	$(CC) $(CFLAGS) -c $< -o $@

# Cleanup
.PHONY: clean
clean:
	rm -rf $(OBJDIR) $(TARGET) $(TEST_TARGET)