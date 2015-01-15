#-----------------------------------------------------------------------
# File    : makefile
# Contents: build item and transaction management (on Unix systems)
# Author  : Christian Borgelt
# History : 2008.10.05 file created from apriori makefile
#           2009.10.09 module pfxtree added
#           2010.02.10 module pfxtree replaced by module clomax
#           2010.10.08 changed standard from -ansi to -std=c99
#           2011.05.06 changed to double support reporting/recording
#           2011.08.24 main program fim16 added (mainly for testing)
#           2012.07.25 module tract with write functions added (trawr)
#           2013.02.03 module $(UTILDIR)/tabwrite.o added
#           2013.03.20 extended the requested warnings in CFBASE
#           2013.09.27 module patspec added (pattern spectrum)
#           2014.10.24 some modules compiled also for double support
#-----------------------------------------------------------------------
SHELL   = /bin/bash
THISDIR = ../../tract/src
UTILDIR = ../../util/src
MATHDIR = ../../math/src

CC      = gcc -std=c99
# CC      = g++
CFBASE  = -Wall -Wextra -Wno-unused-parameter -Wconversion \
          -pedantic $(ADDFLAGS)
CFLAGS  = $(CFBASE) -DNDEBUG -O3 -funroll-loops
# CFLAGS  = $(CFBASE) -DNDEBUG -O3 -funroll-loops -DALIGN8
# CFLAGS  = $(CFBASE) -g
# CFLAGS  = $(CFBASE) -g -DSTORAGE
# CFLAGS  = $(CFBASE) -g -DALIGN8
INCS    = -I$(THISDIR) -I$(UTILDIR) -I$(MATHDIR)

LD      = gcc
# LD      = g++
LDFLAGS = $(ADDFLAGS)
LIBS    = -lm

# ADDOBJS = $(UTILDIR)/storage.o

HDRS    = $(UTILDIR)/arrays.h   $(UTILDIR)/memsys.h   \
          $(UTILDIR)/symtab.h   $(UTILDIR)/escape.h   \
          $(UTILDIR)/tabread.h  $(UTILDIR)/tabwrite.h \
          $(UTILDIR)/scanner.h  tract.h patspec.h clomax.h report.h
OBJS    = $(UTILDIR)/arrays.o   $(UTILDIR)/memsys.o   \
          $(UTILDIR)/idmap.o    $(UTILDIR)/escape.o   \
          $(UTILDIR)/tabread.o  $(UTILDIR)/tabwrite.o \
          $(UTILDIR)/scform.o   patspec.o clomax.o repcm.o $(ADDOBJS)
PRGS    = fim16 tract train

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:          $(PRGS)

fim16:        $(OBJS) tract.o m16main.o makefile
	$(LD) $(LDFLAGS) $(OBJS) tract.o m16main.o $(LIBS) -o $@

tract:        $(OBJS) tramain.o makefile
	$(LD) $(LDFLAGS) $(OBJS) tramain.o $(LIBS) -o $@

train:        $(OBJS) $(UTILDIR)/random.o tract.o trnmain.o makefile
	$(LD) $(LDFLAGS) $(OBJS) $(UTILDIR)/random.o tract.o \
              trnmain.o $(LIBS) -o $@

psp:          $(UTILDIR)/arrays.o   $(UTILDIR)/escape.o   \
              $(UTILDIR)/idmap.o    $(UTILDIR)/tabread.o  \
              $(UTILDIR)/tabwrite.o tract.o train.o
psp:          pspmain.o makefile
	$(LD) $(LDFLAGS) \
              $(UTILDIR)/arrays.o   $(UTILDIR)/escape.o   \
              $(UTILDIR)/idmap.o    $(UTILDIR)/tabread.o  \
              $(UTILDIR)/tabwrite.o tract.o train.o       \
              pspmain.o $(LIBS) -o $@

rgt:          $(UTILDIR)/arrays.o   $(UTILDIR)/escape.o   \
              $(UTILDIR)/idmap.o    $(UTILDIR)/tabread.o  \
              $(UTILDIR)/memsys.o   $(UTILDIR)/scform.o   \
              $(MATHDIR)/ruleval.o  $(MATHDIR)/gamma.o    \
              $(MATHDIR)/chi2.o     \
              tract.o report.o patspec.o
rgt:          rgmain.o makefile
	$(LD) $(LDFLAGS) \
              $(UTILDIR)/arrays.o   $(UTILDIR)/escape.o   \
              $(UTILDIR)/idmap.o    $(UTILDIR)/tabread.o  \
              $(UTILDIR)/memsys.o   $(UTILDIR)/scform.o   \
              $(MATHDIR)/ruleval.o  $(MATHDIR)/gamma.o    \
              $(MATHDIR)/chi2.o     \
              tract.o report.o patspec.o rgmain.o $(LIBS) -o $@

#-----------------------------------------------------------------------
# Main Programs
#-----------------------------------------------------------------------
tramain.o:     $(HDRS)
tramain.o:    tract.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTA_READ -DTA_MAIN -c tract.c -o $@

trnmain.o:    $(HDRS) $(UTILDIR)/random.h $(UTILDIR)/tabwrite.h
trnmain.o:    train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTRN_READ -DTRN_MAIN -c train.c -o $@

m16main.o:    $(HDRS)
m16main.o:    fim16.c fim16.h makefile
	$(CC) $(CFLAGS) $(INCS) -DM16_MAIN -c fim16.c -o $@

pspmain.o:    $(HDRS) $(UTILDIR)/tabwrite.h
pspmain.o:    patspec.c makefile
	$(CC) $(CFLAGS) $(INCS) -DPSP_MAIN -c patspec.c -o $@

rgmain.o:     $(HDRS) $(UTILDIR)/memsys.h $(MATHDIR)/ruleval.h
rgmain.o:     rulegen.c makefile
	$(CC) $(CFLAGS) $(INCS) -DRG_MAIN -c rulegen.c -o $@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tract.o:      tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                      $(UTILDIR)/tabread.h
tract.o:      tract.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTA_READ -c tract.c -o $@

tads.o:       tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                      $(UTILDIR)/tabread.h
tads.o:       tract.c makefile
	$(CC) $(CFLAGS) $(INCS) -DSUPP=double -DTA_READ -c tract.c -o $@

tasurr.o:     tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                      $(UTILDIR)/tabread.h
tasurr.o:     tract.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTA_SURR -c tract.c -o $@

trawr.o:      tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                      $(UTILDIR)/tabread.h $(UTILDIR)/tabwrite.h
trawr.o:      tract.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTA_READ -DTA_WRITE -c tract.c -o $@

tatree.o:     tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h
tatree.o:     tract.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTA_READ -DTATREEFN -c tract.c -o $@

#-----------------------------------------------------------------------
# Train Management
#-----------------------------------------------------------------------
train.o:      train.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h \
                              $(UTILDIR)/tabread.h
train.o:      train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTRN_READ -c train.c -o $@

trnds.o:      train.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h \
                              $(UTILDIR)/tabread.h
trnds.o:      train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DSUPP=double -DTRN_READ \
              -c train.c -o $@

trnwr.o:      train.h tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                              $(UTILDIR)/tabread.h $(UTILDIR)/tabwrite.h
trnwr.o:      train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTRN_READ -DTRN_WRITE -c train.c -o $@

trnsur.o:     train.h tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                              $(UTILDIR)/tabread.h $(UTILDIR)/tabwrite.h
trnsur.o:     train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTRN_SURR -c train.c -o $@

trnrs.o:      train.h tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                              $(UTILDIR)/tabread.h $(UTILDIR)/tabwrite.h
trnrs.o:      train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTRN_READ -DTRN_SURR -c train.c -o $@

trnrws.o:     train.h tract.h $(UTILDIR)/arrays.h  $(UTILDIR)/symtab.h \
                              $(UTILDIR)/tabread.h $(UTILDIR)/tabwrite.h
trnrws.o:     train.c makefile
	$(CC) $(CFLAGS) $(INCS) -DTRN_READ -DTRN_WRITE -DTRN_SURR \
              -c train.c -o $@

#-----------------------------------------------------------------------
# Frequent Item Set Mining (with at most 16 items)
#-----------------------------------------------------------------------
fim16.o:      $(HDRS)
fim16.o:      fim16.c makefile
	$(CC) $(CFLAGS) $(INCS) -c fim16.c -o $@

#-----------------------------------------------------------------------
# Pattern Statistics Management
#-----------------------------------------------------------------------
patspec.o:    patspec.h $(UTILDIR)/tabwrite.h
patspec.o:    patspec.c makefile
	$(CC) $(CFLAGS) $(INCS) -DPSP_REPORT -c patspec.c -o $@

pspdbl.o:     patspec.h $(UTILDIR)/tabwrite.h
pspdbl.o:     patspec.c makefile
	$(CC) $(CFLAGS) $(INCS) -DSUPP=double -DPSP_REPORT \
              -c patspec.c -o $@

pspest.o:     patspec.h $(UTILDIR)/tabwrite.h
pspest.o:     patspec.c makefile
	$(CC) $(CFLAGS) $(INCS) -DPSP_REPORT -DPSP_ESTIM -DPSP_TRAIN \
              -c patspec.c -o $@

#-----------------------------------------------------------------------
# Closed and Maximal Frequent Item Set Tree Management
#-----------------------------------------------------------------------
clomax.o:     clomax.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/memsys.h
clomax.o:     clomax.c makefile
	$(CC) $(CFLAGS) $(INCS) -c clomax.c -o $@

cmdbl.o:      clomax.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/memsys.h
cmdbl.o:      clomax.c makefile
	$(CC) $(CFLAGS) $(INCS) -DRSUPP=double -c clomax.c -o $@

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
report.o:     report.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h
report.o:     report.c makefile
	$(CC) $(CFLAGS) $(INCS) -DISR_PATSPEC -c report.c -o $@

repdbl.o:     report.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h
repdbl.o:     report.c makefile
	$(CC) $(CFLAGS) $(INCS) -DRSUPP=double -DISR_PATSPEC \
              -c report.c -o $@

repcm.o:      report.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h
repcm.o:      report.c makefile
	$(CC) $(CFLAGS) $(INCS) -DISR_PATSPEC -DISR_CLOMAX \
              -c report.c -o $@

repcmd.o:     report.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/symtab.h
repcmd.o:     report.c makefile
	$(CC) $(CFLAGS) $(INCS) -DRSUPP=double -DISR_PATSPEC \
              -DISR_CLOMAX -c report.c -o $@

#-----------------------------------------------------------------------
# Rule Generation Tree Management
#-----------------------------------------------------------------------
rulegen.o:    rulegen.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/memsys.h
rulegen.o:    rulegen.c makefile
	$(CC) $(CFLAGS) $(INCS) -c rulegen.c -o $@

rgrfn.o:      rulegen.h tract.h $(UTILDIR)/arrays.h $(UTILDIR)/memsys.h
rgrfn.o:      rulegen.c makefile
	$(CC) $(CFLAGS) $(INCS) -DRG_REPOFN -c rulegen.c -o $@

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
$(UTILDIR)/storage.o:
	cd $(UTILDIR);  $(MAKE) storage.o  ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/ruleval.o:
	cd $(MATHDIR);  $(MAKE) ruleval.o  ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/gamma.o:
	cd $(MATHDIR);  $(MAKE) gamma.o    ADDFLAGS="$(ADDFLAGS)"
$(MATHDIR)/chi2.o:
	cd $(MATHDIR);  $(MAKE) chi2.o     ADDFLAGS="$(ADDFLAGS)"

#-----------------------------------------------------------------------
# Source Distribution Packages
#-----------------------------------------------------------------------
dist:
	$(MAKE) clean
	cd ../..; rm -f eclat.zip eclat.tar.gz; \
        zip -rq eclat.zip eclat/{src,ex,doc} \
          tract/src/{tract.[ch],train.[ch],report.[ch],clomax.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],makefile,math.mak} math/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],random.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],scanner.[ch],random.[ch]} \
          util/src/{makefile,util.mak} util/doc; \
        tar cfz eclat.tar.gz eclat/{src,ex,doc} \
          tract/src/{tract.[ch],train.[ch],report.[ch],clomax.[ch]} \
          tract/src/{makefile,tract.mak} tract/doc \
          math/src/{gamma.[ch],makefile,math.mak} math/doc \
          util/src/{fntypes.h,error.h} \
          util/src/{arrays.[ch],random.[ch],memsys.[ch],symtab.[ch]} \
          util/src/{escape.[ch],tabread.[ch],scanner.[ch],random.[ch]} \
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
	rm -f *.o *~ *.flc core $(PRGS) psp rgt

clean:
	$(MAKE) localclean
	cd $(UTILDIR); $(MAKE) clean
	cd $(MATHDIR); $(MAKE) clean
