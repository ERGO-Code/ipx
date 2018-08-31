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
	$(CXX) $(CF) $(SO_OPTS) $^ -o $@ $(LDLIBS)
	( cd lib; ln -sf $(SO_TARGET) $(SO_PLAIN) )
	( cd lib; ln -sf $(SO_TARGET) $(SO_MAIN) )

#-------------------------------------------------------------------------------
# compile each object file from its source file
#-------------------------------------------------------------------------------

build/%.o: src/%.cc $(DEP_FILES)
	$(CXX) $(CF) $(IFLAGS) -c $< -o $@

#-------------------------------------------------------------------------------
# compile examples
#-------------------------------------------------------------------------------

examples:
	@( cd example; $(MAKE); )

#-------------------------------------------------------------------------------
# clean and purge
#-------------------------------------------------------------------------------

clean:
	$(RM) $(OBJ_FILES)
	@( cd example; $(MAKE) clean; )

purge: clean
	$(RM) lib/$(AR_TARGET) lib/$(SO_TARGET) lib/$(SO_PLAIN) lib/$(SO_MAIN)
	@( cd example; $(MAKE) purge; )
