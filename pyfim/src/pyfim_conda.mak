#-----------------------------------------------------------------------
# File    : pyfim_conda.mak
# Contents: build pyfim dynamic link library (on Windows systems)
#           (for Anaconda distribution)
# Author  : Christian Borgelt
# History : 2012.08.03 file created
#           2013.10.19 module patspec added
#           2013.11.08 properly adapted to Windows/Microsoft Visual C
#           2013.11.14 adapted for Anaconda Python distribution
#-----------------------------------------------------------------------
!IFNDEF CONDAINC
CONDAINC = "C:\Anaconda\include"
!ENDIF
!IFNDEF CONDALIB
CONDALIB = "C:\Anaconda\libs"
!ENDIF
THISDIR  = ..\..\pyfim\src
UTILDIR  = ..\..\util\src
MATHDIR  = ..\..\math\src
TRACTDIR = ..\..\tract\src
APRIDIR  = ..\..\apriori\src
ECLATDIR = ..\..\eclat\src
FPGDIR   = ..\..\fpgrowth\src
SAMDIR   = ..\..\sam\src
RELIMDIR = ..\..\relim\src
CARPDIR  = ..\..\carpenter\src
ACCDIR   = ..\..\accretion\src

CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_WARNINGS
CFLAGS   = /nologo /MD /W3 /O2 /GS- $(DEFS) /c
PYINC    = /I $(CONDAINC)
INCS     = /I $(UTILDIR) /I $(MATHDIR)  /I $(TRACTDIR) \
           /I $(APRIDIR) /I $(ECLATDIR) /I $(FPGDIR)   \
		   /I $(SAMDIR)  /I $(RELIMDIR) /I $(CARPDIR)  \
           /I $(ACCDIR)

LD       = link.exe
LDFLAGS  = /DLL /nologo /incremental:no
LIBS     = /LIBPATH:$(CONDALIB)

HDRS     = $(UTILDIR)\arrays.h    $(UTILDIR)\memsys.h  \
           $(UTILDIR)\symtab.h    $(MATHDIR)\gamma.h   \
           $(MATHDIR)\chi2.h      $(MATHDIR)\ruleval.h \
           $(TRACTDIR)\tract.h    $(TRACTDIR)\fim16.h  \
           $(TRACTDIR)\patspec.h  $(TRACTDIR)\clomax.h \
           $(TRACTDIR)\report.h
OBJS     = arrays.obj memsys.obj idmap.obj \
           chi2.obj gamma.obj ruleval.obj tatree.obj \
           fim16.obj patspec.obj clomax.obj report.obj \
           istree.obj apriori.obj eclat.obj fpgrowth.obj \
           sam.obj relim.obj repotree.obj carpenter.obj \
		   accretion.obj pyfim.obj
PRGS     = fim.pyd

#-----------------------------------------------------------------------
# Build Dynamic Link Library
#-----------------------------------------------------------------------
all:          $(PRGS)

fim.pyd:      $(OBJS) pyfim.mak
	$(LD) $(LDFLAGS) $(LIBS) $(OBJS) /out:$@ /IMPLIB:fim.lib

#-----------------------------------------------------------------------
# Array Operations
#-----------------------------------------------------------------------
arrays.obj:   $(UTILDIR)\arrays.h $(UTILDIR)\fntypes.h
arrays.obj:   $(UTILDIR)\arrays.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)\arrays.c /Fo$@

#-----------------------------------------------------------------------
# Memory Management System for Objects of Equal Size
#-----------------------------------------------------------------------
memsys.obj:   $(UTILDIR)\memsys.h
memsys.obj:   $(UTILDIR)\memsys.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(UTILDIR)\memsys.c /Fo$@

#-----------------------------------------------------------------------
# Symbol Table Management
#-----------------------------------------------------------------------
idmap.obj:    $(UTILDIR)\symtab.h $(UTILDIR)\fntypes.h \
              $(UTILDIR)\arrays.h
idmap.obj:    $(UTILDIR)\symtab.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D IDMAPFN $(UTILDIR)\symtab.c /Fo$@

#-----------------------------------------------------------------------
# Mathematical Functions
#-----------------------------------------------------------------------
gamma.obj:    $(MATHDIR)\gamma.h
gamma.obj:    $(MATHDIR)\gamma.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)\gamma.c /Fo$@

chi2.obj:     $(MATHDIR)\chi2.h
chi2.obj:     $(MATHDIR)\chi2.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)\chi2.c /Fo$@

ruleval.obj:  $(MATHDIR)\ruleval.h
ruleval.obj:  $(MATHDIR)\ruleval.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(MATHDIR)\ruleval.c /Fo$@

#-----------------------------------------------------------------------
# Item and Transaction Management
#-----------------------------------------------------------------------
tatree.obj:   $(TRACTDIR)\tract.h $(UTILDIR)\arrays.h \
              $(UTILDIR)\symtab.h
tatree.obj:   $(TRACTDIR)\tract.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D TATREEFN $(TRACTDIR)\tract.c /Fo$@

#-----------------------------------------------------------------------
# Item Set Reporter Management
#-----------------------------------------------------------------------
patspec.obj:  $(TRACTDIR)\patspec.h $(TRACTDIR)\tract.h
patspec.obj:  $(TRACTDIR)\patspec.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)\patspec.c /Fo$@

clomax.obj:   $(TRACTDIR)\clomax.h $(TRACTDIR)\tract.h \
              $(UTILDIR)\arrays.h
clomax.obj:   $(TRACTDIR)\clomax.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)\clomax.c /Fo$@

report.obj:   $(TRACTDIR)\report.h  $(TRACTDIR)\clomax.h \
              $(TRACTDIR)\patspec.h $(TRACTDIR)\tract.h \
              $(UTILDIR)\arrays.h   $(UTILDIR)\symtab.h
report.obj:   $(TRACTDIR)\report.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) /D ISR_PATSPEC /D ISR_CLOMAX \
              /D ISR_NONAMES $(TRACTDIR)\report.c /Fo$@

#-----------------------------------------------------------------------
# 16 Items Machine
#-----------------------------------------------------------------------
fim16.obj:    $(TRACTDIR)\tract.h $(TRACTDIR)\report.h \
              $(UTILDIR)\arrays.h  $(UTILDIR)\symtab.h
fim16.obj:    $(TRACTDIR)\fim16.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(TRACTDIR)\fim16.c /Fo$@

#-----------------------------------------------------------------------
# Accretion
#-----------------------------------------------------------------------
accretion.obj: $(HDRS) $(ACCDIR)\accretion.h $(UTILDIR)\fntypes.h
accretion.obj: $(ACCDIR)\accretion.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(ACCDIR)\accretion.c /Fo$@

#-----------------------------------------------------------------------
# Apriori
#-----------------------------------------------------------------------
istree.obj:   $(HDRS) $(APRIDIR)\istree.h
istree.obj:   $(APRIDIR)\istree.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(APRIDIR)\istree.c /Fo$@

apriori.obj:  $(HDRS) $(APRIDIR)\istree.h $(APRIDIR)\apriori.h \
              $(UTILDIR)\fntypes.h
apriori.obj:  $(APRIDIR)\apriori.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(APRIDIR)\apriori.c /Fo$@

#-----------------------------------------------------------------------
# Eclat
#-----------------------------------------------------------------------
eclat.obj:    $(HDRS) $(ECLATDIR)\eclat.h $(UTILDIR)\fntypes.h
eclat.obj:    $(ECLATDIR)\eclat.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(ECLATDIR)\eclat.c /Fo$@

#-----------------------------------------------------------------------
# FP-growth
#-----------------------------------------------------------------------
fpgrowth.obj: $(HDRS) $(FPGDIR)\fpgrowth.h $(UTILDIR)\fntypes.h
fpgrowth.obj: $(FPGDIR)\fpgrowth.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(FPGDIR)\fpgrowth.c /Fo$@

#-----------------------------------------------------------------------
# SaM
#-----------------------------------------------------------------------
sam.obj:      $(HDRS) $(SAMDIR)\sam.h $(UTILDIR)\fntypes.h
sam.obj:      $(SAMDIR)\sam.c makefile
	$(CC) $(CFLAGS) $(INCS) $(SAMDIR)\sam.c /Fo$@

#-----------------------------------------------------------------------
# RElim
#-----------------------------------------------------------------------
relim.obj:    $(HDRS) $(RELIMDIR)\relim.h $(UTILDIR)\fntypes.h
relim.obj:    $(RELIMDIR)\relim.c makefile
	$(CC) $(CFLAGS) $(INCS) $(RELIMDIR)\relim.c /Fo$@

#-----------------------------------------------------------------------
# Carpenter
#-----------------------------------------------------------------------
repotree.obj: $(HDRS) $(CARPDIR)\repotree.h
repotree.obj: $(CARPDIR)\repotree.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(CARPDIR)\repotree.c /Fo$@

carpenter.obj: $(HDRS) $(CARPDIR)\carpenter.h $(CARPDIR)\repotree.h \
               $(UTILDIR)\fntypes.h
carpenter.obj: $(CARPDIR)\carpenter.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(CARPDIR)\carpenter.c /Fo$@

#-----------------------------------------------------------------------
# Python Stuff
#-----------------------------------------------------------------------
pyfim.obj:      $(HDRS) $(ACCDIR)\accretion.h $(APRIDIR)\apriori.h \
                $(ECLATDIR)\eclat.h $(FPGDIR)\fpgrowth.h \
                $(CARPDIR)\carpenter.h
pyfim.obj:      pyfim.c pyfim.mak
	$(CC) $(CFLAGS) $(INCS) $(PYINC) pyfim.c /Fo$@

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	-@copy $(PRGS) ..\..\..\bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
clean:
	-@erase /Q *~ *.obj *.idb *.pch *.exp *.lib $(PRGS)
