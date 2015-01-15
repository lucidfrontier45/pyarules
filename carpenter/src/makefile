#-----------------------------------------------------------------------
# File    : makefile
# Contents: build carpenter program (on Unix systems)
# Author  : Christian Borgelt
# History : 2010.06.23 file created from eclat makefile
#           2010.07.22 repotree module added (specialized clomax)
#           2010.08.22 module escape added (for module tabread)
#           2010.10.08 changed standard from -ansi to -std=c99
#           2013.03.20 extended the requested warnings in CFBASE
#           2013.10.19 modules tabread and patspec added
#-----------------------------------------------------------------------
# For large file support (> 2GB) compile with
#   make ADDFLAGS=-D_FILE_OFFSET_BITS=64
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../carpenter/src
UTILDIR  = ../../util/src
TRACTDIR = ../../tract/src

CC       = gcc -std=c99
# CC       = g++
CFBASE   = -Wall -Wextra -Wno-unused-parameter -Wconversion \
           -pedantic $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
INCS     = -I$(UTILDIR) -I$(TRACTDIR)

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/arrays.h   $(UTILDIR)/memsys.h   \
           $(UTILDIR)/symtab.h   $(UTILDIR)/escape.h   \
           $(UTILDIR)/tabread.h  $(UTILDIR)/tabwrite.h \
           $(UTILDIR)/scanner.h  $(TRACTDIR)/tract.h   \
           $(TRACTDIR)/patspec.h $(TRACTDIR)/clomax.h  \
           $(TRACTDIR)/report.h  repotree.h 
OBJS     = $(UTILDIR)/arrays.o   $(UTILDIR)/memsys.o   \
           $(UTILDIR)/idmap.o    $(UTILDIR)/escape.o   \
           $(UTILDIR)/tabread.o  $(UTILDIR)/tabwrite.o \
           $(UTILDIR)/scform.o   $(TRACTDIR)/tract.o   \
           $(TRACTDIR)/patspec.o $(TRACTDIR)/clomax.o \
           $(TRACTDIR)/repcm.o   repotree.o \
           carpenter.o $(ADDOBJS)
PRGS     = carpenter

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:       $(PRGS)

carpenter: $(OBJS) makefile
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
carpenter.o:  $(HDRS)
carpenter.o:  carpenter.c makefile
	$(CC) $(CFLAGS) $(INCS) -D CARP_MAIN -c carpenter.c -o $@

#-----------------------------------------------------------------------
# Item Set Repository Tree Management
#-----------------------------------------------------------------------
repotree.o:   repotree.h $(TRACTDIR)/report.h $(UTILDIR)/memsys.h
repotree.o:   repotree.c makefile
	$(CC) $(CFLAGS) $(INCS) -c repotree.c -o $@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)/arrays.o:
	cd $(UTILDIR);  $(MAKE) arrays.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/memsys.o:
	cd $(UTILDIR);  $(MAKE) memsys.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/idmap.o:
	cd $(UTILDIR);  $(MAKE) idmap.o    ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/escape.o:
	cd $(UTILDIR);  $(MAKE) escape.o   ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabread.o:
	cd $(UTILDIR);  $(MAKE) tabread.o  ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/tabwrite.o:
	cd $(UTILDIR);  $(MAKE) tabwrite.o ADDFLAGS="$(ADDFLAGS)"
$(UTILDIR)/scform.o:
	cd $(UTILDIR);  $(MAKE) scform.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/tract.o:
	cd $(TRACTDIR); $(MAKE) tract.o    ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/patspec.o:
	cd $(TRACTDIR); $(MAKE) patspec.o  ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/clomax.o:
	cd $(TRACTDIR); $(MAKE) clomax.o   ADDFLAGS="$(ADDFLAGS)"
$(TRACTDIR)/repcm.o:
	cd $(TRACTDIR); $(MAKE) repcm.o    ADDFLAGS="$(ADDFLAGS)"

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f carpenter.zip carpenter.tar.gz; \
        zip -rq carpenter.zip carpenter/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz carpenter.tar.gz carpenter/{src,ex,doc} \
          tract/src/{tract.[ch],patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],tabwrite.[ch]} \
          util/src/scanner.[ch] \
          util/src/{makefile,util.mak} util/doc

#-----------------------------------------------------------------------
# Installation
#-----------------------------------------------------------------------
install:
	cp $(PRGS) $(HOME)/bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	rm -f *.o *~ *.flc core $(PRGS)

clean:
	$(MAKE) localclean
	cd $(TRACTDIR); $(MAKE) localclean
	cd $(UTILDIR);  $(MAKE) clean
