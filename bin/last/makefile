CXXFLAGS = -O3 -std=c++11 -pthread -DHAS_CXX_THREADS
all:
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp src/last?? src/last-split src/last-merge-batches src/last-pair-probs scripts/* $(bindir)

clean:
	@cd src && $(MAKE) clean

html:
	@cd doc && $(MAKE)

distdir = last-`hg id -n`

RSYNCFLAGS = -aC --exclude 'last??' --exclude last-split --exclude last-merge-batches --exclude last-pair-probs

dist: log html
	@cd src && $(MAKE) version.hh CyclicSubsetSeedData.hh ScoreMatrixData.hh
	rsync $(RSYNCFLAGS) build doc examples makefile scripts src data *.txt $(distdir)
	zip -qrm $(distdir) $(distdir)

log:
	hg log --style changelog > ChangeLog.txt
