#===============================================================================
# ipx build configuration
#===============================================================================

# The file has been copied and adapted from SuiteSparse_config.mk contained
# in SuiteSparse version 4.5.3 and available from http://www.suitesparse.com.
# As with the original file, no licensing restrictions apply to this file.

IPX_VERSION = 1.0.0

# Path to the top-level directory of the basiclu package, which can be obtained
# from https://github.com/ERGO-Code/basiclu. By default assume that ipx/ and
# basiclu/ reside side-by-side. If basiclu/include and basiclu/lib are in the
# compiler's default search path, then the line can be commented out.
BASICLUROOT = $(realpath ../basiclu)

#===============================================================================
# Defaults for any system
#===============================================================================

    #---------------------------------------------------------------------------
    # ipx root directory
    #---------------------------------------------------------------------------

    IPXROOT = $(realpath $(CURDIR))

    #---------------------------------------------------------------------------
    # optimization level
    #---------------------------------------------------------------------------

    OPTIMIZATION ?= -O2

    #---------------------------------------------------------------------------
    # compiler flags for both C and C++ compiler
    #---------------------------------------------------------------------------

    COMPILEFLAGS = $(CPPFLAGS) $(TARGET_ARCH) $(OPTIMIZATION) -fPIC

    #---------------------------------------------------------------------------
    # C compiler, must support C99
    #---------------------------------------------------------------------------

    CC = gcc
    CCFLAGS = -std=c99

    #---------------------------------------------------------------------------
    # C++ compiler, must support C++11
    #---------------------------------------------------------------------------

    CXX = g++
    CXXFLAGS = -std=c++11 -Wall -Wextra -Wno-sign-compare

    #---------------------------------------------------------------------------
    # required libraries
    #---------------------------------------------------------------------------

    LAPACK ?= -llapack
    BLAS ?= -lopenblas

    # uncomment if your BLAS library uses 64 bit integers rather than 'int'
    # COMPILEFLAGS += -DBLAS64

    LDLIBS =
    LDLIBS += -lbasiclu
    LDLIBS += $(LAPACK) $(BLAS)

    #---------------------------------------------------------------------------
    # include flags and linker options
    #---------------------------------------------------------------------------

    COMPILEFLAGS += -Iinclude -Isrc
    LINKFLAGS = -Wl,-Bstatic -lipx -lbasiclu -Llib
    SO_OPTS =

    ifdef BASICLUROOT
        COMPILEFLAGS += -I$(BASICLUROOT)/include
        LINKFLAGS += -L$(BASICLUROOT)/lib
        SO_OPTS += -L$(BASICLUROOT)/lib -Wl,-rpath,$(BASICLUROOT)/lib
    endif
    # If BASICLUROOT is not defined, assume that it is in the default search
    # paths.

    LINKFLAGS += -Wl,-Bdynamic
    LINKFLAGS += $(LAPACK) $(BLAS)
    LINKFLAGS += -lstdc++ -lm

    #---------------------------------------------------------------------------
    # shell commands
    #---------------------------------------------------------------------------

    # ranlib, and ar, for generating libraries.  If you don't need ranlib,
    # just change it to RANLIB = echo
    RANLIB = ranlib
    ARCHIVE = $(AR) $(ARFLAGS)

#===============================================================================
# System-dependent configurations
#===============================================================================

    #---------------------------------------------------------------------------
    # determine what system we are on
    #---------------------------------------------------------------------------

    # To disable these auto configurations, use 'make UNAME=custom'

    ifndef UNAME
        ifeq ($(OS),Windows_NT)
            # Cygwin Make on Windows has an $(OS) variable, but not uname.
            # Note that this option is untested.
            UNAME = Windows
        else
            # Linux and Darwin (Mac OSX) have been tested.
            UNAME := $(shell uname)
        endif
    endif

#===============================================================================
# Building the shared and static libraries
#===============================================================================

LIBRARY = libipx
VERSION = 1.0.0
SO_VERSION = 1

ifeq ($(UNAME),Windows)
    # Cygwin / Mingw Make on Windows
    AR_TARGET = $(LIBRARY).lib
    SO_PLAIN  = $(LIBRARY).dll
    SO_MAIN   = $(LIBRARY).$(SO_VERSION).dll
    SO_TARGET = $(LIBRARY).$(VERSION).dll
    SO_OPTS   += -shared
    SO_INSTALL_NAME = echo
else
    # Mac or Linux/Unix
    AR_TARGET = $(LIBRARY).a
    ifeq ($(UNAME),Darwin)
        # Mac
        SO_PLAIN  = $(LIBRARY).dylib
        SO_MAIN   = $(LIBRARY).$(SO_VERSION).dylib
        SO_TARGET = $(LIBRARY).$(VERSION).dylib
        SO_OPTS   += -dynamiclib -compatibility_version $(SO_VERSION) \
                     -current_version $(VERSION) \
                     -shared -undefined dynamic_lookup
        # When a Mac *.dylib file is moved, this command is required
        # to change its internal name to match its location in the filesystem:
        SO_INSTALL_NAME = install_name_tool -id
    else
        # Linux and other variants of Unix
        SO_PLAIN  = $(LIBRARY).so
        SO_MAIN   = $(LIBRARY).so.$(SO_VERSION)
        SO_TARGET = $(LIBRARY).so.$(VERSION)
        SO_OPTS   += -shared -Wl,-soname -Wl,$(SO_MAIN) -Wl,--no-undefined
        # Linux/Unix *.so files can be moved without modification:
        SO_INSTALL_NAME = echo
    endif
endif
