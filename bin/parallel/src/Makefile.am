bin_SCRIPTS = parallel sql niceload parcat env_parallel		\
	env_parallel.ash env_parallel.bash env_parallel.csh	\
	env_parallel.dash env_parallel.fish env_parallel.ksh	\
	env_parallel.pdksh env_parallel.sh env_parallel.tcsh	\
	env_parallel.zsh

install-exec-hook:
	rm $(DESTDIR)$(bindir)/sem || true
	$(LN_S) parallel $(DESTDIR)$(bindir)/sem

if DOCUMENTATION
man_MANS = parallel.1 env_parallel.1 sem.1 sql.1 niceload.1 \
	parallel_tutorial.7 parallel_design.7 parallel_alternatives.7 \
	parcat.1 parset.1
doc_DATA = parallel.html env_parallel.html sem.html sql.html niceload.html \
	parallel_tutorial.html parallel_design.html parallel_alternatives.html \
	 parcat.html \
	parallel.texi env_parallel.texi sem.texi sql.texi niceload.texi \
	parallel_tutorial.texi parallel_design.texi parallel_alternatives.texi \
	parcat.texi parset.texi \
	parallel.pdf env_parallel.pdf sem.pdf sql.pdf niceload.pdf \
	parallel_tutorial.pdf parallel_design.pdf parallel_alternatives.pdf \
	parcat.pdf parset.pdf
endif

parset.pod: parset
	cp parset parset.pod

# Build documentation file if the tool to build exists.
# Otherwise: Use the distributed version
parallel.1: parallel.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/parallel.pod > $(srcdir)/parallel.1n \
	&& mv $(srcdir)/parallel.1n $(srcdir)/parallel.1 \
	|| echo "Warning: pod2man not found. Using old parallel.1"

env_parallel.1: env_parallel.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/env_parallel.pod > $(srcdir)/env_parallel.1n \
	&& mv $(srcdir)/env_parallel.1n $(srcdir)/env_parallel.1 \
	|| echo "Warning: pod2man not found. Using old env_parallel.1"

parallel_tutorial.7: parallel_tutorial.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=7 $(srcdir)/parallel_tutorial.pod > $(srcdir)/parallel_tutorial.7n \
	&& mv $(srcdir)/parallel_tutorial.7n $(srcdir)/parallel_tutorial.7 \
	|| echo "Warning: pod2man not found. Using old parallel_tutorial.7"

parallel_design.7: parallel_design.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=7 $(srcdir)/parallel_design.pod > $(srcdir)/parallel_design.7n \
	&& mv $(srcdir)/parallel_design.7n $(srcdir)/parallel_design.7 \
	|| echo "Warning: pod2man not found. Using old parallel_design.7"


parallel_alternatives.7: parallel_alternatives.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=7 $(srcdir)/parallel_alternatives.pod > $(srcdir)/parallel_alternatives.7n \
	&& mv $(srcdir)/parallel_alternatives.7n $(srcdir)/parallel_alternatives.7 \
	|| echo "Warning: pod2man not found. Using old parallel_alternatives.7"

sem.1: sem.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/sem.pod > $(srcdir)/sem.1n \
	&& mv $(srcdir)/sem.1n $(srcdir)/sem.1 \
	|| echo "Warning: pod2man not found. Using old sem.1"

sql.1: sql
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/sql > $(srcdir)/sql.1n \
	&& mv $(srcdir)/sql.1n $(srcdir)/sql.1 \
	|| echo "Warning: pod2man not found. Using old sql.1"

niceload.1: niceload.pod
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/niceload.pod > $(srcdir)/niceload.1n \
	&& mv $(srcdir)/niceload.1n $(srcdir)/niceload.1 \
	|| echo "Warning: pod2man not found. Using old niceload.1"

parcat.1: parcat
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/parcat > $(srcdir)/parcat.1n \
	&& mv $(srcdir)/parcat.1n $(srcdir)/parcat.1 \
	|| echo "Warning: pod2man not found. Using old parcat.1"

parset.1: parset
	pod2man --release='$(PACKAGE_VERSION)' --center='$(PACKAGE_NAME)' \
		--section=1 $(srcdir)/parset > $(srcdir)/parset.1n \
	&& mv $(srcdir)/parset.1n $(srcdir)/parset.1 \
	|| echo "Warning: pod2man not found. Using old parset.1"

parallel.html: parallel.pod
	pod2html --title "GNU Parallel" $(srcdir)/parallel.pod > $(srcdir)/parallel.htmln \
	&& mv $(srcdir)/parallel.htmln $(srcdir)/parallel.html \
	|| echo "Warning: pod2html not found. Using old parallel.html"
	rm -f $(srcdir)/pod2htm*

# Depending on parallel.html to avoid stupid pod2html race condition
env_parallel.html: env_parallel.pod parallel.html
	pod2html --title "GNU Parallel with environment" $(srcdir)/env_parallel.pod > $(srcdir)/env_parallel.htmln \
	&& mv $(srcdir)/env_parallel.htmln $(srcdir)/env_parallel.html \
	|| echo "Warning: pod2html not found. Using old env_parallel.html"
	rm -f $(srcdir)/pod2htm*

# Depending on env_parallel.html to avoid stupid pod2html race condition
parallel_tutorial.html: parallel_tutorial.pod env_parallel.html
	pod2html --title "GNU Parallel tutorial" $(srcdir)/parallel_tutorial.pod > $(srcdir)/parallel_tutorial.htmln \
	&& mv $(srcdir)/parallel_tutorial.htmln $(srcdir)/parallel_tutorial.html \
	|| echo "Warning: pod2html not found. Using old parallel_tutorial.html"
	rm -f $(srcdir)/pod2htm*

# Depending on parallel_tutorial.html to avoid stupid pod2html race condition
parallel_design.html: parallel_design.pod parallel_tutorial.html
	pod2html --title "GNU Parallel design" $(srcdir)/parallel_design.pod > $(srcdir)/parallel_design.htmln \
	&& mv $(srcdir)/parallel_design.htmln $(srcdir)/parallel_design.html \
	|| echo "Warning: pod2html not found. Using old parallel_design.html"
	rm -f $(srcdir)/pod2htm*

# Depending on parallel_design.html to avoid stupid pod2html race condition
parallel_alternatives.html: parallel_alternatives.pod parallel_design.html
	pod2html --title "GNU Parallel alternatives" $(srcdir)/parallel_alternatives.pod > $(srcdir)/parallel_alternatives.htmln \
	&& mv $(srcdir)/parallel_alternatives.htmln $(srcdir)/parallel_alternatives.html \
	|| echo "Warning: pod2html not found. Using old parallel_alternatives.html"
	rm -f $(srcdir)/pod2htm*

# Depending on parallel_alternatives.html to avoid stupid pod2html race condition
sem.html: sem.pod parallel_alternatives.html
	pod2html --title "sem (GNU Parallel)" $(srcdir)/sem.pod > $(srcdir)/sem.htmln \
	&& mv $(srcdir)/sem.htmln $(srcdir)/sem.html \
	|| echo "Warning: pod2html not found. Using old sem.html"
	rm -f $(srcdir)/pod2htm*

# Depending on sem.html to avoid stupid pod2html race condition
sql.html: sql sem.html
	pod2html --title "GNU SQL" $(srcdir)/sql > $(srcdir)/sql.htmln \
	&& mv $(srcdir)/sql.htmln $(srcdir)/sql.html \
	|| echo "Warning: pod2html not found. Using old sql.html"
	rm -f $(srcdir)/pod2htm*

# Depending on sql.html to avoid stupid pod2html race condition
niceload.html: niceload.pod sql.html
	pod2html --title "GNU niceload" $(srcdir)/niceload.pod > $(srcdir)/niceload.htmln \
	&& mv $(srcdir)/niceload.htmln $(srcdir)/niceload.html \
	|| echo "Warning: pod2html not found. Using old niceload.html"
	rm -f $(srcdir)/pod2htm*

# Depending on niceload.html to avoid stupid pod2html race condition
parcat.html: parcat niceload.html
	pod2html --title "GNU parcat" $(srcdir)/parcat > $(srcdir)/parcat.htmln \
	&& mv $(srcdir)/parcat.htmln $(srcdir)/parcat.html \
	|| echo "Warning: pod2html not found. Using old parcat.html"
	rm -f $(srcdir)/pod2htm*

# Depending on niceload.html to avoid stupid pod2html race condition
parset.html: parset parcat.html
	pod2html --title "GNU parset" $(srcdir)/parset > $(srcdir)/parset.htmln \
	&& mv $(srcdir)/parset.htmln $(srcdir)/parset.html \
	|| echo "Warning: pod2html not found. Using old parset.html"
	rm -f $(srcdir)/pod2htm*

parallel.texi: parallel.pod
	pod2texi --output=$(srcdir)/parallel.texi $(srcdir)/parallel.pod \
	|| echo "Warning: pod2texi not found. Using old parallel.texi"

env_parallel.texi: env_parallel.pod
	pod2texi --output=$(srcdir)/env_parallel.texi $(srcdir)/env_parallel.pod \
	|| echo "Warning: pod2texi not found. Using old env_parallel.texi"

parallel_tutorial.texi: parallel_tutorial.pod
	pod2texi --output=$(srcdir)/parallel_tutorial.texi $(srcdir)/parallel_tutorial.pod \
	|| echo "Warning: pod2texi not found. Using old parallel_tutorial.texi"

parallel_design.texi: parallel_design.pod
	pod2texi --output=$(srcdir)/parallel_design.texi $(srcdir)/parallel_design.pod \
	|| echo "Warning: pod2texi not found. Using old parallel_design.texi"

parallel_alternatives.texi: parallel_alternatives.pod
	pod2texi --output=$(srcdir)/parallel_alternatives.texi $(srcdir)/parallel_alternatives.pod \
	|| echo "Warning: pod2texi not found. Using old parallel_alternatives.texi"

sem.texi: sem.pod
	pod2texi --output=$(srcdir)/sem.texi $(srcdir)/sem.pod \
	|| echo "Warning: pod2texi not found. Using old sem.texi"

sql.texi: sql
	pod2texi --output=$(srcdir)/sql.texi $(srcdir)/sql \
	|| echo "Warning: pod2texi not found. Using old sql.texi"

niceload.texi: niceload.pod
	pod2texi --output=$(srcdir)/niceload.texi $(srcdir)/niceload.pod \
	|| echo "Warning: pod2texi not found. Using old niceload.texi"

parcat.texi: parcat
	pod2texi --output=$(srcdir)/parcat.texi $(srcdir)/parcat \
	|| echo "Warning: pod2texi not found. Using old parcat.texi"

parset.texi: parset
	pod2texi --output=$(srcdir)/parset.texi $(srcdir)/parset \
	|| echo "Warning: pod2texi not found. Using old parset.texi"

parallel.pdf: parallel.pod
	pod2pdf --output-file $(srcdir)/parallel.pdf $(srcdir)/parallel.pod --title "GNU Parallel" \
	|| echo "Warning: pod2pdf not found. Using old parallel.pdf"

env_parallel.pdf: env_parallel.pod
	pod2pdf --output-file $(srcdir)/env_parallel.pdf $(srcdir)/env_parallel.pod --title "GNU Parallel with environment" \
	|| echo "Warning: pod2pdf not found. Using old env_parallel.pdf"

parallel_tutorial.pdf: parallel_tutorial.pod
	pod2pdf --output-file $(srcdir)/parallel_tutorial.pdf $(srcdir)/parallel_tutorial.pod --title "GNU Parallel Tutorial" \
	|| echo "Warning: pod2pdf not found. Using old parallel_tutorial.pdf"

parallel_design.pdf: parallel_design.pod
	pod2pdf --output-file $(srcdir)/parallel_design.pdf $(srcdir)/parallel_design.pod --title "GNU Parallel Design" \
	|| echo "Warning: pod2pdf not found. Using old parallel_design.pdf"

parallel_alternatives.pdf: parallel_alternatives.pod
	pod2pdf --output-file $(srcdir)/parallel_alternatives.pdf $(srcdir)/parallel_alternatives.pod --title "GNU Parallel alternatives" \
	|| echo "Warning: pod2pdf not found. Using old parallel_alternatives.pdf"

sem.pdf: sem.pod
	pod2pdf --output-file $(srcdir)/sem.pdf $(srcdir)/sem.pod --title "GNU sem" \
	|| echo "Warning: pod2pdf not found. Using old sem.pdf"

sql.pdf: sql
	pod2pdf --output-file $(srcdir)/sql.pdf $(srcdir)/sql --title "GNU SQL" \
	|| echo "Warning: pod2pdf not found. Using old sql.pdf"

niceload.pdf: niceload.pod
	pod2pdf --output-file $(srcdir)/niceload.pdf $(srcdir)/niceload.pod --title "GNU niceload" \
	|| echo "Warning: pod2pdf not found. Using old niceload.pdf"

parcat.pdf: parcat
	pod2pdf --output-file $(srcdir)/parcat.pdf $(srcdir)/parcat --title "GNU parcat" \
	|| echo "Warning: pod2pdf not found. Using old parcat.pdf"

parset.pdf: parset
	pod2pdf --output-file $(srcdir)/parset.pdf $(srcdir)/parset --title "GNU parset" \
	|| echo "Warning: pod2pdf not found. Using old parset.pdf"

sem: parallel
	ln -fs parallel sem

DISTCLEANFILES = parallel.1 env_parallel.1 sem.1 sql.1 niceload.1 \
	parallel_tutorial.7 parallel_design.7 parallel_alternatives.7 \
	parcat.1 parset.1 \
	parallel.html env_parallel.html sem.html sql.html niceload.html \
	parallel_tutorial.html parallel_design.html parallel_alternatives.html \
	parcat.html parset.html \
	parallel.texi env_parallel.texi sem.texi sql.texi niceload.texi \
	parallel_tutorial.texi parallel_design.texi parallel_alternatives.texi \
	parcat.texi parset.texi \
	parallel.pdf env_parallel.pdf sem.pdf sql.pdf niceload.pdf \
	parallel_tutorial.pdf parallel_design.pdf parallel_alternatives.pdf \
	parcat.pdf parset.pdf

EXTRA_DIST = parallel sem sql niceload parcat parset env_parallel	\
	env_parallel.ash env_parallel.bash env_parallel.csh		\
	env_parallel.dash env_parallel.fish env_parallel.ksh		\
	env_parallel.pdksh env_parallel.sh env_parallel.tcsh		\
	env_parallel.zsh sem.pod parallel.pod env_parallel.pod		\
	niceload.pod parallel_tutorial.pod parallel_design.pod		\
	parallel_alternatives.pod $(DISTCLEANFILES)
