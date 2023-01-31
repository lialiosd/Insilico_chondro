CXX      := icx
CXXFLAGS := -Wall -L$(MKLROOT)/lib/intel64
CPPFLAGS := -Iinclude -I$(MKLROOT)/include
OPTRPT   :=  #-qopt-report=5
##LDLIBS   := -L$(LD_LIBRARY_PATH) -L$(MKLROOT)/lib/intel64
LDFLAGS  := -lstdc++ -qopenmp -qmkl=parallel

SRC_DIR := ./src
OBJ_DIR := ./obj
SRC := $(wildcard $(SRC_DIR)/*.cc)
OBJ := $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

default : app

# Link the executable

.PHONY :app
app : $(OBJ)
	$(CXX) $(LDFLAGS) $(LDLIBS) $^ -o $@

# Build specific files

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc
	$(CXX) $(CPPFLAGS) $(CXX_FLAGS) -c $< -o $@

.PHONY :clean
clean :
	rm app *.o *.optrpt ./obj/*
