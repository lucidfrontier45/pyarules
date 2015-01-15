#-----------------------------------------------------------------------
# File    : ista.mak
# Contents: build ista program (intersect transactions)(Windows systems)
# Author  : Christian Borgelt
# History : 2009.10.15 file created from sam makefile
#           2010.08.22 module escape added (for module tabread)
#           2011.12.01 module clomax added (for module report)
#-----------------------------------------------------------------------
THISDIR  = ..\..\ista\src
UTILDIR  = ..\..\util\src
TRACTDIR = ..\..\tract\src

CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_WARNINGS
CFLAGS   = /nologo /W3 /O2 /GS- $(DEFS) /c $(ADDFLAGS)
INCS     = /I $(UTILDIR) /I $(TRACTDIR)

LD       = link.exe
LDFLAGS  = /nologo /subsystem:console /incremental:no
LIBS     = 

HDRS     = $(UTILDIR)\arrays.h     $(UTILDIR)\memsys.h     \
           $(UTILDIR)\symtab.h     $(UTILDIR)\escape.h     \
           $(UTILDIR)\tabread.h    $(UTILDIR)\tabwrite.h   \
           $(UTILDIR)\scanner.h    $(TRACTDIR)\tract.h     \
           $(TRACTDIR)\patspec.h   $(TRACTDIR)\clomax.h    \
           $(TRACTDIR)\report.h    pfxtree.h pattree.h
OBJS     = $(UTILDIR)\arrays.obj   $(UTILDIR)\memsys.obj   \
           $(UTILDIR)\idmap.obj    $(UTILDIR)\escape.obj   \
           $(UTILDIR)\tabread.obj  $(UTILDIR)\tabwrite.obj \
           $(UTILDIR)\scform.obj   $(TRACTDIR)\tract.obj   \
           $(TRACTDIR)\patspec.obj $(TRACTDIR)\clomax.obj  \
           $(TRACTDIR)\repcm.obj   pfxtree.obj pattree.obj ista.obj
PRGS     = ista.exe

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:       $(PRGS)

ista.exe:  $(OBJS) ista.mak
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) /out:$@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
ista.obj:  $(HDRS) ista.mak
	$(CC) $(CFLAGS) $(INCS) /D ISTA_MAIN ista.c /Fo$@

#-----------------------------------------------------------------------
# Prefix Tree Management
#-----------------------------------------------------------------------
pfxtree.obj: $(UTILDIR)\memsys.h pfxtree.h pfxtree.c ista.mak
	$(CC) $(CFLAGS) $(INCS) pfxtree.c /Fo$@

#-----------------------------------------------------------------------
# Patricia Tree Management
#-----------------------------------------------------------------------
pattree.obj: $(UTILDIR)\memsys.h pattree.h pattree.c ista.mak
	$(CC) $(CFLAGS) $(INCS) pattree.c /Fo$@

#-----------------------------------------------------------------------
# External Modules
#-----------------------------------------------------------------------
$(UTILDIR)\arrays.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak arrays.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\memsys.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak memsys.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\idmap.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak idmap.obj    ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\escape.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak escape.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\tabread.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak tabread.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\tabwrite.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak tabwrite.obj ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(UTILDIR)\scform.obj:
	cd $(UTILDIR)
	$(MAKE) /f util.mak scform.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\tract.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak tract.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\patspec.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak patspec.obj ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\clomax.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak clomax.obj  ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)
$(TRACTDIR)\repcm.obj:
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak repcm.obj   ADDFLAGS="$(ADDFLAGS)"
	cd $(THISDIR)

#-----------------------------------------------------------------------
# Install
#-----------------------------------------------------------------------
install:
	-@copy $(PRGS) ..\..\..\bin

#-----------------------------------------------------------------------
# Clean up
#-----------------------------------------------------------------------
localclean:
	-@erase /Q *~ *.obj *.idb *.pch $(PRGS)

clean:
	$(MAKE) /f ista.mak localclean
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak localclean
	cd $(UTILDIR)
	$(MAKE) /f util.mak clean
	cd $(THISDIR)
