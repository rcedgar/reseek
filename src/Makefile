BINDIR := ../bin
OBJDIR := o
BINPATH := $(BINDIR)/reseek

CXX = ccache g++
CXXFLAGS := $(CFLAGS) -DNDEBUG -pthread -O3 -ffast-math -march=native --std=c++11

UNAME_S := $(shell uname -s)
LDFLAGS := $(LDFLAGS) -O3 -pthread -lpthread
ifeq ($(UNAME_S),Linux)
    LDFLAGS += -static
endif

HDRS = \
  abcxyz.h \
  active.h \
  alpha.h \
  best.h \
  binner.h \
  calibratesearcher.h \
  calreader.h \
  chainreader.h \
  cigar.h \
  cigar_template.h \
  cmds.h \
  combomx.h \
  counters.h \
  dbsearcher.h \
  derepper.h \
  dss.h \
  dssaligner.h \
  dssp.h \
  dssparams.h \
  fast_erf.h \
  fastaseqsource.h \
  featurelist.h \
  features.h \
  fileseqsource.h \
  filler.h \
  floatfeatures.h \
  getticks.h \
  gobuff.h \
  intfeatures.h \
  linereader.h \
  lockobj.h \
  lockobjs.h \
  logodds.h \
  megacluster.h \
  murmur.h \
  mx.h \
  myopts.h \
  myutils.h \
  obj.h \
  objmgr.h \
  objtype.h \
  objtypes.h \
  omplock.h \
  parasail.h \
  pathinfo.h \
  pdbchain.h \
  quarts.h \
  scalarparams.h \
  scop40bench.h \
  seqdb.h \
  seqinfo.h \
  seqsource.h \
  sfasta.h \
  sort.h \
  spher.h \
  sweeper.h \
  swtrace.h \
  timers.h \
  timing.h \
  tma.h \
  tracebit.h \
  trainer.h \
  usage.h \
  xdpmem.h \

OBJS = \
  $(OBJDIR)/abcxyz.o \
  $(OBJDIR)/alpha.o \
  $(OBJDIR)/alpha2.o \
  $(OBJDIR)/b62.o \
  $(OBJDIR)/blosum.o \
  $(OBJDIR)/blosum62.o \
  $(OBJDIR)/calibrate.o \
  $(OBJDIR)/calreader.o \
  $(OBJDIR)/chainreader.o \
  $(OBJDIR)/cigar.o \
  $(OBJDIR)/explore1.o \
  $(OBJDIR)/explore1i.o \
  $(OBJDIR)/parasail.o \
  $(OBJDIR)/parasail_cigar.o \
  $(OBJDIR)/parasail_combo.o \
  $(OBJDIR)/pdb2mega.o \
  $(OBJDIR)/readhits.o \
  $(OBJDIR)/scop40bit.o \
  $(OBJDIR)/search.o \
  $(OBJDIR)/combomx.o \
  $(OBJDIR)/combomx_data.o \
  $(OBJDIR)/dali.o \
  $(OBJDIR)/dbsearcher.o \
  $(OBJDIR)/derepper.o \
  $(OBJDIR)/dss.o \
  $(OBJDIR)/dssaligner.o \
  $(OBJDIR)/dssparams.o \
  $(OBJDIR)/eval_feature.o \
  $(OBJDIR)/featurefreqs.o \
  $(OBJDIR)/features.o \
  $(OBJDIR)/float_feature_bins.o \
  $(OBJDIR)/getbins.o \
  $(OBJDIR)/gsw.o \
  $(OBJDIR)/jiggle.o \
  $(OBJDIR)/kabsch.o \
  $(OBJDIR)/kmer_eval.o \
  $(OBJDIR)/kmer_eval2.o \
  $(OBJDIR)/logaln.o \
  $(OBJDIR)/logodds.o \
  $(OBJDIR)/mergefwdback.o \
  $(OBJDIR)/msa2cmp.o \
  $(OBJDIR)/murmurhash.o \
  $(OBJDIR)/myss.o \
  $(OBJDIR)/namedparams.o \
  $(OBJDIR)/pdb2cal.o \
  $(OBJDIR)/prettyaln.o \
  $(OBJDIR)/scalar_sweep.o \
  $(OBJDIR)/scop40bench.o \
  $(OBJDIR)/scop40benchroc.o \
  $(OBJDIR)/scoremxs2.o \
  $(OBJDIR)/sweeper.o \
  $(OBJDIR)/swfastpinop.o \
  $(OBJDIR)/swfastpinopgapless.o \
  $(OBJDIR)/swgapless.o \
  $(OBJDIR)/swgaplessint.o \
  $(OBJDIR)/swgaplessprof.o \
  $(OBJDIR)/swgaplessprofb.o \
  $(OBJDIR)/test_para.o \
  $(OBJDIR)/test_para_cal.o \
  $(OBJDIR)/test_para_path.o \
  $(OBJDIR)/test_xdrop.o \
  $(OBJDIR)/timing.o \
  $(OBJDIR)/usort.o \
  $(OBJDIR)/valuetoint.o \
  $(OBJDIR)/one2three.o \
  $(OBJDIR)/pdb2dss.o \
  $(OBJDIR)/fastaseqsource.o \
  $(OBJDIR)/fileseqsource.o \
  $(OBJDIR)/getss.o \
  $(OBJDIR)/sscluster.o \
  $(OBJDIR)/trainer.o \
  $(OBJDIR)/usage.o \
  $(OBJDIR)/linereader.o \
  $(OBJDIR)/lockobj.o \
  $(OBJDIR)/mx.o \
  $(OBJDIR)/myutils.o \
  $(OBJDIR)/objmgr.o \
  $(OBJDIR)/reseek_main.o \
  $(OBJDIR)/pdbchaincal.o \
  $(OBJDIR)/readchains.o \
  $(OBJDIR)/pdbchain.o \
  $(OBJDIR)/quarts.o \
  $(OBJDIR)/spher.o \
  $(OBJDIR)/sw.o \
  $(OBJDIR)/seqdb.o \
  $(OBJDIR)/seqinfo.o \
  $(OBJDIR)/seqsource.o \
  $(OBJDIR)/sfasta.o \
  $(OBJDIR)/sort.o \
  $(OBJDIR)/xdropbwd.o \
  $(OBJDIR)/xdropfwd.o \

.PHONY: clean

$(BINPATH) : $(BINDIR)/ $(OBJDIR)/ $(OBJS)
	$(CXX) $(LDFLAGS) $(OBJS) -o $(BINPATH)
	strip -d $(BINPATH)

$(OBJDIR)/ :
	mkdir -p $(OBJDIR)/

$(BINDIR)/ :
	mkdir -p $(BINDIR)/

$(OBJDIR)/%.o : %.cpp $(HDRS)
	$(CXX) $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<

clean:
	rm -rf $(OBJDIR)/ $(BINPATH)
