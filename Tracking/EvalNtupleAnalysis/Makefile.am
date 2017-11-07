AUTOMAKE_OPTIONS = foreign

INCLUDES = \
  -I$(includedir) \
  -I$(OFFLINE_MAIN)/include \
  -I$(ROOTSYS)/include \
  -I`root-config --incdir`

lib_LTLIBRARIES = libevalanal.la

noinst_HEADERS = \
  EvalNtupleAnalysis.h	EvalNtupleAnalysisLinkDef.h 

include_HEADERS = \
  EvalNtupleAnalysis.h  

libevalanal_la_LDFLAGS = \
  -nodefaultlibs \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib \
  `root-config --libs`

libevalanal_la_SOURCES = \
  EvalNtupleAnalysis.C	EvalNtupleAnalysisDict.C 

libevalanal_la_LIBADD = \
  -L$(libdir) \
  -L$(OFFLINE_MAIN)/lib

# Rule for generating CINT dictionaries from class headers.
%Dict.C: %.h %LinkDef.h
	rootcint -f $@ -c $(DEFAULT_INCLUDES) $(INCLUDES) $^

clean-local:
	rm -f *Dict*

