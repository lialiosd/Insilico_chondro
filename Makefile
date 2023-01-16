CXX      := icpc
CXXFLAGS := -Wall #-qopenmp
CPPFLAGS := -Iinclude
OPTRPT   :=  #-qopt-report=5
LDLIBS   := -lm
LDFLAGS  := -Llib

SRC_DIR := ./src
OBJ_DIR := ./obj
SRC := $(wildcard $(SRC_DIR)/*.cc)
OBJ := $(SRC:$(SRC_DIR)/%.cc=$(OBJ_DIR)/%.o)

default : app


# Build the executable

.PHONY :app
app : $(OBJ)
	$(CXX) $(LDFLAGS) $^ $(LDLIBS) -o $@

# Build specific files

$(OBJ_DIR)/%.o : $(SRC_DIR)/%.cc
	$(CXX) $(CPPFLAGS) $(CXX_FLAGS) -c $< -o $@

.PHONY :clean
clean :
	rm app *.o *.optrpt
