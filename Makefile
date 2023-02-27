include config.mk

#-------------------------------------------------------------------------------
# compile the package
#-------------------------------------------------------------------------------

# The default creates the static and shared library in lib/ and compiles the
# *.c and *.cc files from example/ into executables stored in bin/example/.

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
EXAMPLE_BIN_FILES := $(patsubst %, bin/%, $(EXAMPLE_BIN_FILES))

TEST_SRC_FILES = $(wildcard check/*.cc)
TEST_DEP_FILES = $(wildcard check/*.h)
TEST_BIN_FILE = bin/$(TEST_TARGET)

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
	$(CXX) $(CXXFLAGS) $(COMPILEFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# compile examples
#-------------------------------------------------------------------------------

examples: $(EXAMPLE_BIN_FILES)

bin/example/%: example/%.cc static
	$(CXX) $(CXXFLAGS) $(COMPILEFLAGS) $< -o $@ $(LINKFLAGS)

bin/example/%: example/%.c static
	$(CC) $(CCFLAGS) $(COMPILEFLAGS) $< -o $@ $(LINKFLAGS)

#-------------------------------------------------------------------------------
# compile and run tests
#-------------------------------------------------------------------------------

test: $(TEST_BIN_FILE)

$(TEST_BIN_FILE): $(TEST_SRC_FILES) $(TEST_DEP_FILES) static
	$(CXX) $(CXXFLAGS) $(COMPILEFLAGS) $(TEST_SRC_FILES) -o $@ $(LINKFLAGS)
	./$@

#-------------------------------------------------------------------------------
# clean and purge
#-------------------------------------------------------------------------------

clean:
	$(RM) $(OBJ_FILES)

purge: clean
	$(RM) lib/$(AR_TARGET) lib/$(SO_TARGET) lib/$(SO_PLAIN) lib/$(SO_MAIN)
	$(RM) $(EXAMPLE_BIN_FILES)
	$(RM) $(TEST_BIN_FILE)
