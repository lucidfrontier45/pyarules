#-----------------------------------------------------------------------
# File    : carpenter.mak
# Contents: build carpenter program (on Windows systems)
# Author  : Christian Borgelt
# History : 2010.06.23 file created from eclat makefile
#           2010.08.22 module escape added (for module tabread)
#           2011.12.01 module clomax added (for module report)
#           2013.10.19 modules tabread and patspec added
#-----------------------------------------------------------------------
THISDIR  = ..\..\carpenter\src
UTILDIR  = ..\..\util\src
TRACTDIR = ..\..\tract\src

CC       = cl.exe
DEFS     = /D WIN32 /D NDEBUG /D _CONSOLE /D _CRT_SECURE_NO_DEPRECATE
CFLAGS   = /nologo /W3 /O2 /GS- $(DEFS) /c $(ADDFLAGS)
INCS     = /I $(UTILDIR) /I $(TRACTDIR)

LD       = link.exe
LDFLAGS  = /nologo /subsystem:console /incremental:no
LIBS     = 

HDRS     = $(UTILDIR)\arrays.h     $(UTILDIR)\memsys.h   \
           $(UTILDIR)\symtab.h     $(UTILDIR)\escape.h   \
           $(UTILDIR)\tabread.h    $(UTILDIR)\tabwrite.h \
           $(UTILDIR)\scanner.h    $(TRACTDIR)\tract.h   \
           $(TRACTDIR)\patspec.h   $(TRACTDIR)\clomax.h  \
           $(TRACTDIR)\report.h    repotree.h
OBJS     = $(UTILDIR)\arrays.obj   $(UTILDIR)\memsys.obj   \
           $(UTILDIR)\idmap.obj    $(UTILDIR)\escape.obj   \
           $(UTILDIR)\tabread.obj  $(UTILDIR)\tabwrite.obj \
           $(UTILDIR)\scform.obj   $(TRACTDIR)\tract.obj   \
           $(TRACTDIR)\patspec.obj $(TRACTDIR)\clomax.obj  \
           $(TRACTDIR)\repcm.obj   repotree.obj carpenter.obj
PRGS     = carpenter.exe

#-----------------------------------------------------------------------
# Build Program
#-----------------------------------------------------------------------
all:         $(PRGS)

carpenter.exe:  $(OBJS) carpenter.mak
	$(LD) $(LDFLAGS) $(OBJS) $(LIBS) /out:$@

#-----------------------------------------------------------------------
# Main Program
#-----------------------------------------------------------------------
carpenter.obj:  $(HDRS) carpenter.mak
	$(CC) $(CFLAGS) $(INCS) /D CARP_MAIN carpenter.c /Fo$@

#-----------------------------------------------------------------------
# Item Set Repository Tree Management
#-----------------------------------------------------------------------
repotree.obj:  $(TRACTDIR)\report.h $(UTILDIR)\memsys.h \
               repotree.h repotree.c carpenter.mak
	$(CC) $(CFLAGS) $(INCS) repotree.c /Fo$@

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
	$(MAKE) /f carpenter.mak localclean
	cd $(TRACTDIR)
	$(MAKE) /f tract.mak localclean
	cd $(UTILDIR)
	$(MAKE) /f util.mak clean
	cd $(THISDIR)
