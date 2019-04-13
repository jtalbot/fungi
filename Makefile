# This Makefile requires GNU make.
 
UNAME := $(shell uname -s)
CXX := clang++
CC := clang
CXXFLAGS := -Wall -std=c++14 `pkg-config --cflags OpenEXR luajit libpng`
CFLAGS := -Wall
LIBS := `pkg-config --libs OpenEXR luajit libpng`
LFLAGS :=

ifeq ($(UNAME),Linux)
	LIBS += -lrt
endif

SRC := fungi.cpp api.cpp pbrt.cpp sensor.cpp material.cpp light.cpp
TEST_SRC := fungi_tests.cpp la_tests.cpp

EXECUTABLE := fungi
TEST_EXECUTABLE := fungi_tests

OBJECTS := $(patsubst %.cpp,build/%.o,$(SRC))
DEPENDENCIES := $(patsubst %.cpp,build/%.d,$(SRC))

TEST_OBJECTS := $(patsubst %.cpp,build/%.o,$(TEST_SRC))
TEST_DEPENDENCIES := $(patsubst %.cpp,build/%.d,$(TEST_SRC))

default: release 

debug: CXXFLAGS += -DDEBUG -O0 -g
debug: all

release: CXXFLAGS += -DNDEBUG -O3 -g
release: all

tests: $(TEST_EXECUTABLE)
	./fungi_tests

all: $(EXECUTABLE) $(TEST_EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LFLAGS) -o $@ $^ $(LIBS)

$(TEST_EXECUTABLE): $(TEST_OBJECTS) 
	$(CXX) $(LFLAGS) -o $@ $^ $(LIBS)

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

build/%.o: tests/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@
 
.PHONY: clean
clean:
	rm -rf $(EXECUTABLE) $(OBJECTS) $(DEPENDENCIES) $(TEST_OBJECTS) $(TEST_DEPENDENCIES)

# dependency rules
build/%.d: src/%.cpp
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT '$@ $(@:.d=.o)' $< -o $@
	
build/%.d: tests/%.cpp
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT '$@ $(@:.d=.o)' $< -o $@
	
-include $(DEPENDENCIES)
-include $(TEST_DEPENDENCIES)

