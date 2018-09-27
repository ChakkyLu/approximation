# simple makefile

MKDIR_P = mkdir -p

# main directories
BIN_DIR = ../bin
OBJ_DIR = ./build
SRC_DIR = ./src
INC_DIR = ./include
LIB_DIR = ../libs

ABC_DIR = ../abc
GLUCOSE_DIR = ../glucose

# lexer options
LEXER = reflex
LEXER_FLAGS =

PARSER_DIR = $(SRC_DIR)/parser
VG_PARSER_LEX = $(PARSER_DIR)/verilog_gate_parser.l
VG_PARSER_CPP = $(SRC_DIR)/parser/verilog_gate_parser.cpp
VG_PARSER_INC = $(INC_DIR)/parser/verilog_gate_parser.h
BL_PARSER_LEX = $(PARSER_DIR)/blif_parser.l
BL_PARSER_CPP = $(SRC_DIR)/parser/blif_parser.cpp
BL_PARSER_INC = $(INC_DIR)/parser/blif_parser.h
VGL_PARSER_LEX = $(PARSER_DIR)/verilog_glib_parser.l
VGL_PARSER_CPP = $(SRC_DIR)/parser/verilog_glib_parser.cpp
VGL_PARSER_INC = $(INC_DIR)/parser/verilog_glib_parser.h


# sub project directories
INC_DIRS = $(shell find $(INC_DIR) -type d)
INC_DIRS += $(GLUCOSE_DIR) $(ABC_DIR)/src

# compiler and linker options
EXE_NAME = cadex

CXX = g++

CXX_FLAGS = -W -Wall -Wextra -s -O3 -static -std=c++17 -DABC_NAMESPACE=ABC -DLIN64
INC_FLAGS = $(addprefix -I,$(INC_DIRS))
LD_FLAGS = -static
LIB_FLAGS = -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -lreflex -lm -ldl -rdynamic -lz

# collect sources ...
SRCS = $(shell find $(SRC_DIR) -name "*.cpp")
OBJS = $(SRCS:%.cpp=$(OBJ_DIR)/%.o)
#DEPS = $(OBJS:.o=.d)
LIBS = $(LIB_DIR)/lib.a $(LIB_DIR)/libabc.a

# rules for c++ files
$(OBJ_DIR)/%.o: %.cpp
	$(MKDIR_P) $(dir $@)
	$(CXX) $(INC_FLAGS) $(CXX_FLAGS) -c $< -o $@

#rules for target
$(BIN_DIR)/$(EXE_NAME): $(OBJS)
	$(MKDIR_P) $(BIN_DIR)
	$(CXX) $(LD_FLAGS) $(OBJS) $(LIBS) $(LIB_FLAGS) -o $@


.PHONY: abc glucose parser superclean clean

abc:
	cd $(ABC_DIR); make ABC_USE_NAMESPACE=ABC ABC_USE_PIC=1 ABC_USE_NO_PTHREADS=1 ABC_USE_NO_READLINE=1 libabc.a
	cp $(ABC_DIR)/libabc.a ${LIB_DIR}

glucose:
	cd $(GLUCOSE_DIR)/simp; make libr
	cp $(GLUCOSE_DIR)/simp/lib.a $(LIB_DIR)

parser:
	$(MKDIR_P) $(PARSER_DIR)
	$(LEXER) -o $(VG_PARSER_CPP) --header-file=$(VG_PARSER_INC) $(VG_PARSER_LEX)
#	sed -i 's/reflex_code_INITIAL/vg_reflex_code_INITIAL/g' $(VG_PARSER_CPP)
	$(LEXER) -o $(BL_PARSER_CPP) --header-file=$(BL_PARSER_INC) $(BL_PARSER_LEX)
#	sed -i 's/reflex_code_INITIAL/bl_reflex_code_INITIAL/g' $(BL_PARSER_CPP)
	$(LEXER) -o $(VGL_PARSER_CPP) --header-file=$(VGL_PARSER_INC) $(VGL_PARSER_LEX)
#	sed -i 's/reflex_code_INITIAL/vgl_reflex_code_INITIAL/g' $(VGL_PARSER_CPP)	

clean:
	rm -rf $(OBJ_DIR)/*
	rm -rf $(BIN_DIR)/*

superclean: clean
#	rm $(PARSER_CPP)
#	rm $(PARSER_INC)
	rm -rf $(LIB_DIR)/*
	cd $(GLUCOSE_DIR)/simp; make clean
	cd $(ABC_DIR); make clean
