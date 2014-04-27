# This Makefile requires GNU make.
 
UNAME := $(shell uname -s)
CXX := clang++
CC := clang
CXXFLAGS := -Wall -std=c++11
CFLAGS := -Wall
LIBS := `pkg-config --libs OpenEXR lua`
LFLAGS :=

ifeq ($(UNAME),Linux)
	LIBS += -lrt
endif

SRC := fungi.cpp

EXECUTABLE := fungi 

OBJECTS := $(patsubst %.cpp,build/%.o,$(SRC))
DEPENDENCIES := $(patsubst %.cpp,build/%.d,$(SRC))

default: release 

debug: CXXFLAGS += -DDEBUG -O0 -g
debug: all

release: CXXFLAGS += -DNDEBUG -O3 -g
release: all

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS) 
	$(CXX) $(LFLAGS) -o $@ $^ $(LIBS)

build/%.o: src/%.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@ 

.PHONY: clean
clean:
	rm -rf $(EXECUTABLE) $(OBJECTS) $(DEPENDENCIES)

# dependency rules
build/%.d: src/%.cpp
	@mkdir -p $(dir $@)
	@$(CXX) $(CXXFLAGS) -MM -MT '$@ $(@:.d=.o)' $< -o $@
	
-include $(DEPENDENCIES)

