include ../config.mk

#-------------------------------------------------------------------------------
# program sources and generated binaries
#-------------------------------------------------------------------------------

SRC_FILES = $(wildcard *.cc)
BIN_FILES = $(patsubst %.cc, %, $(SRC_FILES))

#-------------------------------------------------------------------------------
# compile all examples
#-------------------------------------------------------------------------------

# The default compiles each *.cc file in this directory into a standalone
# executable. If you want to compile a single program, call e.g. 'make afiro',
# which will compile afiro.cc.

default: $(BIN_FILES)

#-------------------------------------------------------------------------------
# compile each source file into its executable
#-------------------------------------------------------------------------------

# The IPX shared library must have been compiled before. We bake rpath into the
# executables so that they will run regardless if libipx.so is in the runtime
# search path. The disadvantage is that the executables cannot be moved to
# another directory and can only be called from within the example/ directory,
# but for example programs that should be OK.

%: %.cc
	$(CXX) $(CF) -I../include -I../src $< -o $@ -lipx -L../lib -Wl,-rpath,../lib

#-------------------------------------------------------------------------------
# clean and purge
#-------------------------------------------------------------------------------

clean:

purge: clean
	$(RM) $(BIN_FILES)
