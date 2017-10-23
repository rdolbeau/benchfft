if DOUBLE
if HAVE_ARMPL
PRG=doit doitr
ARMPL_LIBS=@LARMPL@
endif
else
if HAVE_ARMPL
PRG=doit doitr
ARMPL_LIBS=@LARMPL@
endif
endif

INCLUDES = $(INCLBENCH)
BENCHMARK_USEROPT=''

doit_SOURCES=doit.c
doit_LDADD=$(ARMPL_LIBS) $(LIBBENCH) @FLIBS@

doitr_SOURCES=doitr.c
doitr_LDADD=$(ARMPL_LIBS) $(LIBBENCH) @FLIBS@

include ../Makefile.common
