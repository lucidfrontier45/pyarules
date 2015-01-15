#-----------------------------------------------------------------------
# File    : makefile
# Contents: build shared object for frequent item set mining in Python
# Author  : Christian Borgelt
# History : 2011.07.08 file created
#           2011.07.22 module ruleval added
#           2013.10.19 module patspec added
#           2013.11.07 modules tabread and scanner removed
#           2013.08.08 modules gamma, chi2 and ruleval added
#           2014.08.25 module ista added
#-----------------------------------------------------------------------
SHELL    = /bin/bash
THISDIR  = ../../pyfim/src
UTILDIR  = ../../util/src
MATHDIR  = ../../math/src
TRACTDIR = ../../tract/src
APRIDIR  = ../../apriori/src
ECLATDIR = ../../eclat/src
FPGDIR   = ../../fpgrowth/src
SAMDIR   = ../../sam/src
RELIMDIR = ../../relim/src
CARPDIR  = ../../carpenter/src
ISTADIR  = ../../ista/src
ACCDIR   = ../../accretion/src

CC       = gcc -std=c99
# CC       = g++
PYINC    = -I/usr/include/python2.7
# PYINC    = -I/usr/include/python3.3
INCS     = -I$(UTILDIR) -I$(MATHDIR)  -I$(TRACTDIR) \
           -I$(APRIDIR) -I$(ECLATDIR) -I$(FPGDIR)   \
           -I$(SAMDIR)  -I$(RELIMDIR) -I$(CARPDIR)  \
           -I$(ISTADIR) -I$(ACCDIR)
CFBASE   = -Wall -Wextra -Wconversion -Wno-unused-parameter \
           -Wno-overlength-strings -pedantic -fPIC $(ADDFLAGS)
CFLAGS   = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS   = $(CFBASE) -DNDEBUG -O3 -DALIGN8
# CFLAGS   = $(CFBASE) -g
# CFLAGS   = $(CFBASE) -g -DSTORAGE
# CFLAGS   = $(CFBASE) -g -DALIGN8
XCFLAGS  = $(CFLAGS) -fvisibility=hidden

LD       = gcc
LDFLAGS  = $(ADDFLAGS)
LIBS     = -lm

# ADDOBJS  = $(UTILDIR)/storage.o

HDRS     = $(UTILDIR)/arrays.h   $(UTILDIR)/memsys.h    \
           $(UTILDIR)/symtab.h   $(UTILDIR)/random.h    \
           $(MATHDIR)/gamma.h    $(MATHDIR)/chi2.h      \
           $(MATHDIR)/ruleval.h  $(TRACTDIR)/tract.h    \
           $(TRACTDIR)/fim16.h   $(TRACTDIR)/patspec.h  \
           $(TRACTDIR)/clomax.h  $(TRACTDIR)/report.h   \
           $(APRIDIR)/apriori.h  $(ECLATDIR)/eclat.h    \
           $(FPGDIR)/fpgrowth.h  $(SAMDIR)/sam.h        \
           $(RELIMDIR)/relim.h   $(CARPDIR)/carpenter.h \
           $(ISTADIR)/ista.h     $(ACCDIR)/accretion.h
OBJS     = arrays.o memsys.o idmap.o random.o chi2.o gamma.o ruleval.o \
           tatree.o fim16.o patspec.o clomax.o report.o \
           istree.o apriori.o eclat.o fpgrowth.o sam.o relim.o \
           repotree.o carpenter.o pfxtree.o pattree.o ista.o \
           accretion.o pyfim.o $(ADDOBJS)

#-----------------------------------------------------------------------
# Build Shared Library
#-----------------------------------------------------------------------
all:          fim.so

fim.so:       $(OBJS)
	$(LD) $(LDFLAGS) -shared $(OBJS) -o $@

#-----------------------------------------------------------------------
# Array Operations
#-----------------------------------------------------------------------
arrays.o:     $(UTILDIR)/arrays.h $(UTILDIR)/fntypes.h
arrays.o:     $(UTILDIR)/arrays.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(UTILDIR)/arrays.c -o $@

#-----------------------------------------------------------------------
# Memory Management System for Objects of Equal Size
#-----------------------------------------------------------------------
memsys.o:     $(UTILDIR)/memsys.h
memsys.o:     $(UTILDIR)/memsys.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(UTILDIR)/memsys.c -o $@

#-----------------------------------------------------------------------
# Symbol Table Management
#-----------------------------------------------------------------------
idmap.o:      $(UTILDIR)/symtab.h $(UTILDIR)/fntypes.h \
              $(UTILDIR)/arrays.h
idmap.o:      $(UTILDIR)/symtab.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DIDMAPFN -c $(UTILDIR)/symtab.c -o $@

#-----------------------------------------------------------------------
# Random Number Generator Management
#-----------------------------------------------------------------------
random.o:     $(UTILDIR)/random.h
random.o:     $(UTILDIR)/random.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(UTILDIR)/random.c -o $@

#-----------------------------------------------------------------------
# Mathematical Functions
#-----------------------------------------------------------------------
gamma.o:      $(MATHDIR)/gamma.h
gamma.o:      $(MATHDIR)/gamma.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(MATHDIR)/gamma.c -o $@

chi2.o:       $(MATHDIR)/chi2.h
chi2.o:       $(MATHDIR)/chi2.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(MATHDIR)/chi2.c -o $@

ruleval.o:    $(MATHDIR)/ruleval.h
ruleval.o:    $(MATHDIR)/ruleval.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(MATHDIR)/ruleval.c -o $@

#-----------------------------------------------------------------------
# 16 Items Machine
#-----------------------------------------------------------------------
fim16.o:      $(TRACTDIR)/tract.h $(TRACTDIR)/report.h \
              $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h
fim16.o:      $(TRACTDIR)/fim16.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(TRACTDIR)/fim16.c -o $@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tatree.o:     $(TRACTDIR)/tract.h $(UTILDIR)/arrays.h \
              $(UTILDIR)/symtab.h
tatree.o:     $(TRACTDIR)/tract.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DTATREEFN -DTA_SURR \
              -c $(TRACTDIR)/tract.c -o $@

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
patspec.o:    $(TRACTDIR)/patspec.h $(TRACTDIR)/tract.h
patspec.o:    $(TRACTDIR)/patspec.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DPSP_ESTIM \
              -c $(TRACTDIR)/patspec.c -o $@

clomax.o:     $(TRACTDIR)/clomax.h $(TRACTDIR)/tract.h \
              $(UTILDIR)/arrays.h
clomax.o:     $(TRACTDIR)/clomax.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(TRACTDIR)/clomax.c -o $@

report.o:     $(TRACTDIR)/report.h  $(TRACTDIR)/clomax.h \
              $(TRACTDIR)/patspec.h $(TRACTDIR)/tract.h \
              $(UTILDIR)/arrays.h   $(UTILDIR)/symtab.h
report.o:     $(TRACTDIR)/report.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISR_PATSPEC -DISR_CLOMAX \
              -DISR_NONAMES -c $(TRACTDIR)/report.c -o $@

#-----------------------------------------------------------------------
# Accretion
#-----------------------------------------------------------------------
accretion.o:  $(HDRS) $(ACCDIR)/accretion.h $(UTILDIR)/fntypes.h
accretion.o:  $(ACCDIR)/accretion.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISR_CLOMAX \
              -c $(ACCDIR)/accretion.c -o $@

#-----------------------------------------------------------------------
# Apriori
#-----------------------------------------------------------------------
istree.o:     $(HDRS) $(APRIDIR)/istree.h
istree.o:     $(APRIDIR)/istree.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DTATREEFN -c $(APRIDIR)/istree.c -o $@

apriori.o:    $(HDRS) $(APRIDIR)/istree.h $(APRIDIR)/apriori.h \
              $(UTILDIR)/fntypes.h
apriori.o:    $(APRIDIR)/apriori.c makefile
	$(CC) $(XCFLAGS) $(INCS) -DISR_CLOMAX \
              -c $(APRIDIR)/apriori.c -o $@

#-----------------------------------------------------------------------
# Eclat
#-----------------------------------------------------------------------
eclat.o:      $(HDRS) $(ECLATDIR)/eclat.h $(UTILDIR)/fntypes.h
eclat.o:      $(ECLATDIR)/eclat.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(ECLATDIR)/eclat.c -o $@

#-----------------------------------------------------------------------
# FP-growth
#-----------------------------------------------------------------------
fpgrowth.o:   $(HDRS) $(FPGDIR)/fpgrowth.h $(UTILDIR)/fntypes.h
fpgrowth.o:   $(FPGDIR)/fpgrowth.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(FPGDIR)/fpgrowth.c -o $@

#-----------------------------------------------------------------------
# SaM
#-----------------------------------------------------------------------
sam.o:        $(HDRS) $(SAMDIR)/sam.h $(UTILDIR)/fntypes.h
sam.o:        $(SAMDIR)/sam.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(SAMDIR)/sam.c -o $@

#-----------------------------------------------------------------------
# RElim
#-----------------------------------------------------------------------
relim.o:      $(HDRS) $(RELIMDIR)/relim.h $(UTILDIR)/fntypes.h
relim.o:      $(RELIMDIR)/relim.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(RELIMDIR)/relim.c -o $@

#-----------------------------------------------------------------------
# Carpenter
#-----------------------------------------------------------------------
repotree.o:   $(HDRS) $(CARPDIR)/repotree.h
repotree.o:   $(CARPDIR)/repotree.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(CARPDIR)/repotree.c -o $@

carpenter.o:  $(HDRS) $(CARPDIR)/carpenter.h $(CARPDIR)/repotree.h \
              $(UTILDIR)/fntypes.h
carpenter.o:  $(CARPDIR)/carpenter.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(CARPDIR)/carpenter.c -o $@

#-----------------------------------------------------------------------
# IsTa
#-----------------------------------------------------------------------
pfxtree.o:    $(HDRS) $(ISTADIR)/pfxtree.h $(UTILDIR)/memsys.h
pfxtree.o:    $(ISTADIR)/pfxtree.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(ISTADIR)/pfxtree.c -o $@

pattree.o:    $(HDRS) $(ISTADIR)/pattree.h
pattree.o:    $(ISTADIR)/pattree.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(ISTADIR)/pattree.c -o $@

ista.o:       $(HDRS) $(ISTADIR)/pattree.h $(ISTADIR)/pfxtree.h \
              $(ISTADIR)/ista.h
ista.o:       $(ISTADIR)/ista.c makefile
	$(CC) $(XCFLAGS) $(INCS) -c $(ISTADIR)/ista.c -o $@

#-----------------------------------------------------------------------
# Python Stuff
#-----------------------------------------------------------------------
pyfim.o:      $(HDRS) $(ACCDIR)/accretion.h $(APRIDIR)/apriori.h \
              $(ECLATDIR)/eclat.h $(FPGDIR)/fpgrowth.h \
              $(SAMDIR)/sam.h $(RELIMDIR)/relim.h \
              $(CARPDIR)/carpenter.h $(ISTADIR)/ista.h
pyfim.o:      pyfim.c makefile
	$(CC) $(CFLAGS) $(INCS) $(PYINC) -fno-strict-aliasing \
              -c pyfim.c -o $@

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f pyfim.zip pyfim.tar.gz; \
	cp setup_fim.py setup.py; \
        zip -rq pyfim.zip pyfim/{src,ex,doc} \
          apriori/{src,ex,doc} eclat/{src,ex,doc} \
          fpgrowth/{src,ex,doc} sam/{src,ex,doc} relim/{src,ex,doc} \
          carpenter/{src,ex,doc} ista/{src,ex,doc} \
          accretion/{src,ex,doc} \
          tract/src/{tract.[ch],fim16.[ch]} \
          tract/src/{patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,arrays.[ch],symtab.[ch],memsys.[ch]} \
          util/src/{random.[ch],makefile,util.mak} util/doc \
          setup_fim.py setup.py vs*.bat; \
        tar cfz pyfim.tar.gz pyfim/{src,ex,doc} \
          apriori/{src,ex,doc} eclat/{src,ex,doc} \
          fpgrowth/{src,ex,doc} sam/{src,ex,doc} relim/{src,ex,doc} \
          carpenter/{src,ex,doc} ista/{src,ex,doc} \
          accretion/{src,ex,doc} \
          tract/src/{tract.[ch],fim16.[ch]} \
          tract/src/{patspec.[ch],clomax.[ch],report.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],chi2.[ch],ruleval.[ch]} \
          math/src/{makefile,math.mak} math/doc \
          util/src/{fntypes.h,arrays.[ch],symtab.[ch],memsys.[ch]} \
          util/src/{random.[ch],makefile,util.mak} util/doc \
          setup_fim.py setup.py vs*.bat; \
	rm -f setup.py

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	mv fim.so $(HOME)/lib

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	rm -f *.o *.so *~ *.flc core
