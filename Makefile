#Makefile
#EMTV3D Simulator Makefile
#Author Jianwen Chen  <jwchen@ee.ucla.edu>

NAME=  emtv3d 

include config.mak

BINDIR= bin
INCDIR= src
SRCDIR= src
OBJDIR= obj
LIBDIR= lib
PAPIDIR= /mnt/jc5/CS259/papi

CFLAGS+= -g -pg -ffloat-store
LDFLAGS+= -L$(PAPIDIR) -lpapi -lutil_papi
INCFLAGS= -I$(PAPIDIR) -I$(INCDIR)

SUFFIX=

OBJSUF= .o$(SUFFIX)

SRC=    $(wildcard $(SRCDIR)/*.c)
OBJ=    $(SRC:$(SRCDIR)/%.c=$(OBJDIR)/%.o$(SUFFIX))
BIN=    $(BINDIR)/$(NAME)$(SUFFIX)
LIB=    $(LIBDIR)/libemtv3d.a

default: depend bin tags

dependencies:
	@echo "" >dependencies

clean:
	@echo remove all objects
	@rm -rf $(OBJDIR)/*
	@rm -rf $(OBJDIR)
	@rm -f $(BIN)
	@rm -f ./*~
	@rm -f ./tags

distclean: clean
	rm -f config.mak config.h

tags:
	@echo update tag table
	@ctags -R ./*

bin:    $(OBJ)
	@echo 'creating binary "$(BIN)"'
	@$(CC) $(CFLAGS) $(INCFLAGS) $(OBJ) -o $(BIN) $(LDFLAGS)
	@echo '... done'
	@echo

depend:
	@echo
	@echo 'checking dependencies'
	@echo 'create the obj directory'
	@mkdir -p $(OBJDIR)
	@echo 'create the bin directory'
	@mkdir -p $(BINDIR)
	@echo


$(OBJDIR)/%.o$(SUFFIX): $(SRCDIR)/%.c
	@echo 'compiling object file "$@" ...'
	$(CC) $(CFLAGS) $(INCFLAGS) -c $< -o $@
