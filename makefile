# Variables
CC = g++
VERSION = -std=c++23
OPT = -O3
DEBUG ?= no
MEMCHECK ?=

# Libraries settings
LIB_NAMES := Gillespie Monte_Carlo
LIB_SRC_DIRS := Shared_Objects
LIB_INCLUDE_DIRS := Shared_Objects $(LIB_NAMES)

LIB_SRCSCOMMON := $(foreach dir, $(LIB_SRC_DIRS), $(wildcard $(dir)/*.cpp))
LIB_SRCSGIL := $(foreach dir, $(word 1, $(LIB_NAMES)), $(wildcard $(dir)/*.cpp))
LIB_SRCSMC := $(foreach dir, $(word 2, $(LIB_NAMES)), $(wildcard $(dir)/*.cpp))


LIB_OBJSMC := $(foreach src, $(LIB_SRCSCOMMON) $(LIB_SRCSMC), $(patsubst %.cpp, %.o, $(src)))
LIB_OBJSGIL := $(foreach src, $(LIB_SRCSCOMMON) $(LIB_SRCSGIL), $(patsubst %.cpp, %.o, $(src)))

LIB_HEADERS := $(foreach dir, $(LIB_INCLUDE_DIRS), $(wildcard $(dir)/*.h))
LIB_SO := $(foreach lib, $(LIB_NAMES), $(lib).so)


# Flags
ifeq ($(DEBUG), yes)
    FLAG = -DDEBUG
else
    FLAG =
endif

# Targets
all: $(LIB_SO)

# Compile each cpp file to a .o file
%.o: %.cpp $(LIB_HEADERS)
	$(CC) $(VERSION) -fPIC $(OPT) -c $< -o $@ $(FLAG) $(MEMCHECK)

# Link all the .o files into a shared object file
Gillespie.so:$(LIB_OBJSGIL)
	$(CC) $(OPT) -shared -Wl,-soname,$@ -o $@ $(LIB_OBJSGIL) $(FLAG)
Monte_Carlo.so:$(LIB_OBJSMC)
	$(CC) $(OPT) -shared -Wl,-soname,$@ -o $@ $(LIB_OBJSMC) $(FLAG)

.PHONY: clean EXEC

clean:
	rm -rf $(foreach dir, $(LIB_INCLUDE_DIRS), $(dir)/*.o) *~ $(LIB_SO)

EXEC: $(LIB_OBJS)
	$(foreach lib, $(LIB_NAMES), $(CC) $(OPT) $(LIB_OBJS) $(lib)/$(lib).o -o $(lib).so $(FLAG) $(MEMCHECK);)
