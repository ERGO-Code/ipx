include config.mk

#-------------------------------------------------------------------------------
# compile the package
#-------------------------------------------------------------------------------

# The default creates the static and shared library in the lib/ subdirectory and
# compiles the examples in the example/ subdirectory.

default: static shared examples

#-------------------------------------------------------------------------------
# files
#-------------------------------------------------------------------------------

DEP_FILES = $(wildcard src/*.h) $(wildcard include/*.h) Makefile config.mk
SRC_FILES = $(wildcard src/*.cc)
OBJ_FILES = $(patsubst src/%.cc, build/%.o, $(SRC_FILES))

EXAMPLE_SRC_FILES = $(wildcard example/*.cc) $(wildcard example/*.c)
EXAMPLE_BIN_FILES = $(EXAMPLE_SRC_FILES)
EXAMPLE_BIN_FILES := $(patsubst %.cc, %, $(EXAMPLE_BIN_FILES))
EXAMPLE_BIN_FILES := $(patsubst %.c, %, $(EXAMPLE_BIN_FILES))

#-------------------------------------------------------------------------------
# create the static library
#-------------------------------------------------------------------------------

static: lib/$(AR_TARGET)

lib/$(AR_TARGET): $(OBJ_FILES)
	$(ARCHIVE) $@ $^
	- $(RANLIB) $@

#-------------------------------------------------------------------------------
# create the shared library
#-------------------------------------------------------------------------------

shared: lib/$(SO_TARGET)

lib/$(SO_TARGET): $(OBJ_FILES)
	$(CXX) $(CXXFLAGS) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd lib; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd lib; ln -sf $(SO_TARGET) $(SO_MAIN) )

#-------------------------------------------------------------------------------
# compile each object file from its source file
#-------------------------------------------------------------------------------

build/%.o: src/%.cc $(DEP_FILES)
	$(CXX) $(CXXFLAGS) $(COMPILEFLAGS) $(IFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# compile examples
#-------------------------------------------------------------------------------

examples: $(EXAMPLE_BIN_FILES)

# We bake rpath into the executables so that they will run regardless of whether
# libipx.so is in the runtime search path. The disadvantage is that the
# executables cannot be moved to another directory and can only be called from
# the examples/ or the top-level directory, but for example programs that should
# be OK.

example/%: example/%.cc shared
	$(CXX) $(CXXFLAGS) $(COMPILEFLAGS) -Iinclude -Isrc $< -o $@ -lipx -Llib -Wl,-rpath,lib -Wl,-rpath,../lib

example/%: example/%.c shared
	$(CC) $(CCFLAGS) $(COMPILEFLAGS) -Iinclude $< -o $@ -lipx -Llib -Wl,-rpath,lib -Wl,-rpath,../lib

#-------------------------------------------------------------------------------
# clean and purge
#-------------------------------------------------------------------------------

clean:
	$(RM) $(OBJ_FILES)

purge: clean
	$(RM) lib/$(AR_TARGET) lib/$(SO_TARGET) lib/$(SO_PLAIN) lib/$(SO_MAIN)
	$(RM) $(EXAMPLE_BIN_FILES)
