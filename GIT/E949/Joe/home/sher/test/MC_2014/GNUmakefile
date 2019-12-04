# $Id: GNUmakefile,v 1.2 2003/01/23 15:31:39 maire Exp $
# --------------------------------------------------------------
# GNUmakefile for examples module.  Gabriele Cosmo, 06/04/98.
# --------------------------------------------------------------

name := pienu
G4TARGET := $(name)
G4EXLIB := true

ifndef G4INSTALL
  G4INSTALL = ../../..
endif

.PHONY: all
all: lib bin

#ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += -DG4ANALYSIS_USE_ROOT
#endif

include $(G4INSTALL)/config/architecture.gmk

#ifdef G4ANALYSIS_USE_ROOT
  CPPFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
  LDFLAGS  += $(shell $(ROOTSYS)/bin/root-config --glibs)
  LDFLAGS += -L/home//pienumgr/apps64/xerces-c/lib 
  LDFLAGS += -Wl,-rpath,/home/pienumgr/apps64/xerces-c/lib 
  LDFLAGS += -lxerces-c
  LDFLAGS +=  -Wl,-rpath,/home/pienumgr/apps64/2011/clhep/lib

#endif

ifdef MAKEHBOOK
 CPPFLAGS+=-I$(G4INSTALL)/tests/tools/include -DMAKEHBOOK -DMAKEBOOK
 EXTRALIBS+=-L$(G4INSTALL)/tests/tools/lib/$(G4SYSTEM) -lG4TestTool
 EXTRALIBS+=-L/cern/pro/lib -lpacklib -lnsl
 ifeq ($(G4SYSTEM),Linux-g++)
   EXTRALIBS+=-lg2c -lm -lgcc -lc -lgcc
 endif
endif

include $(G4INSTALL)/config/binmake.gmk

visclean:
	rm -f g4*.prim g4*.eps g4*.wrl
	rm -f .DAWN_*

anaclean:
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/pienuHisto*
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/pienuAnalysis*
	rm -f $(G4WORKDIR)/tmp/$(G4SYSTEM)/$(G4TARGET)/pienuStepping*
