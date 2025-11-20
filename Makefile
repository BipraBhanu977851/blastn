# Makefile for Simple BLASTN Program
# Compiles with C++17 standard and optimization level O2

CXX = g++
CXXFLAGS = -std=c++17 -O2 -Wall
TARGET = simple_blastn
SOURCES = main.cpp fasta.cpp index.cpp search.cpp scoring.cpp
OBJECTS = $(SOURCES:.cpp=.o)

# Default target
all: $(TARGET)

# Link object files to create executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJECTS)

# Compile source files to object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean build artifacts
clean:
	rm -f $(OBJECTS) $(TARGET)

# Rebuild from scratch
rebuild: clean all

# Run the program (builds if needed)
run: $(TARGET)
	./$(TARGET) --db database.fasta --query query.fasta

# Phony targets
.PHONY: all clean rebuild run
