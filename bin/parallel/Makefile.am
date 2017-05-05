SUBDIRS = src

test:
	cd testsuite; make 1

upload:
	eval `gpg-agent --daemon`
# make sure $YYYYMMDD is set
	echo $(YYYYMMDD) | grep 20
	cp parallel-$(YYYYMMDD).tar.bz2 parallel-latest.tar.bz2
	cp doc/parallel.directive parallel-$(YYYYMMDD).tar.bz2.directive
	perl -i -pe "s/20\d\d\d\d\d\d/$(YYYYMMDD)/" parallel-*.tar.*directive
	gpg --clearsign --yes parallel-$(YYYYMMDD).tar.bz2.directive

	cp doc/parallel.latest.directive parallel-latest.tar.bz2.directive
	perl -i -pe "s/20\d\d\d\d\d\d/$(YYYYMMDD)/" parallel-latest.tar.*directive
	gpg --clearsign --yes parallel-latest.tar.bz2.directive

	(echo '#!/bin/bash'; \
	 echo; \
	 echo "# To check the signature run:"; \
	 echo "#  echo | gpg"; \
	 echo "#  gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-$(YYYYMMDD).tar.bz2.sig"; \
	 echo; \
	 echo "echo | gpg 2>/dev/null"; \
	 echo 'gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve $$0'; \
	 echo 'exit $$?' ; \
	 echo; \
	 gpg -ab -o - parallel-$(YYYYMMDD).tar.bz2; \
	) > parallel-$(YYYYMMDD).tar.bz2.sig

	(echo '#!/bin/bash'; \
	 echo; \
	 echo "# To check the signature run:"; \
	 echo "#  echo | gpg"; \
	 echo "#  gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-latest.tar.bz2.sig"; \
	 echo; \
	 echo "echo | gpg 2>/dev/null"; \
	 echo 'gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve $$0'; \
	 echo 'exit $$?' ; \
	 echo; \
	 gpg -ab -o - parallel-$(YYYYMMDD).tar.bz2; \
	) > parallel-latest.tar.bz2.sig

	gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-latest.tar.bz2.sig
	gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-$(YYYYMMDD).tar.bz2.sig

	../ftpsync/src/ftpsync parallel-$(YYYYMMDD).tar.bz2{,.sig,*asc} ftp://ftp-upload.gnu.org/incoming/ftp/
	../ftpsync/src/ftpsync parallel-latest.tar.bz2{,.sig,*asc} ftp://ftp-upload.gnu.org/incoming/ftp/

# This can take 7 minutes
	pushd /tmp; \
	rm -rf /tmp/parallel-$(YYYYMMDD)*; \
	while ! wget http://ftp.gnu.org/gnu/parallel/parallel-$(YYYYMMDD).tar.bz2 ; do sleep 2; done; \
	tar xjvf parallel-$(YYYYMMDD).tar.bz2; \
	cd parallel-$(YYYYMMDD); \
	./configure; \
	make -j && sudo make -j install; \
	pushd; \
	sudo cp /usr/local/bin/parallel /usr/local/bin/parallel-$(YYYYMMDD)

alphaupload:
# Copy of upload
	eval `gpg-agent --daemon`
# make sure $YYYYMMDD is set
	echo $(YYYYMMDD) | grep 20
	cp parallel-$(YYYYMMDD).tar.bz2 parallel-latest.tar.bz2
	cp doc/parallel.directive parallel-$(YYYYMMDD).tar.bz2.directive
	perl -i -pe "s/20\d\d\d\d\d\d/$(YYYYMMDD)/" parallel-*.tar.*directive
	gpg --clearsign --yes parallel-$(YYYYMMDD).tar.bz2.directive

	cp doc/parallel.latest.directive parallel-latest.tar.bz2.directive
	perl -i -pe "s/20\d\d\d\d\d\d/$(YYYYMMDD)/" parallel-latest.tar.*directive
	gpg --clearsign --yes parallel-latest.tar.bz2.directive

	(echo '#!/bin/bash'; \
	 echo; \
	 echo "# To check the signature run:"; \
	 echo "#  echo | gpg"; \
	 echo "#  gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-$(YYYYMMDD).tar.bz2.sig"; \
	 echo; \
	 echo "echo | gpg 2>/dev/null"; \
	 echo 'gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve $$0'; \
	 echo 'exit $$?' ; \
	 echo; \
	 gpg -ab -o - parallel-$(YYYYMMDD).tar.bz2; \
	) > parallel-$(YYYYMMDD).tar.bz2.sig

	(echo '#!/bin/bash'; \
	 echo; \
	 echo "# To check the signature run:"; \
	 echo "#  echo | gpg"; \
	 echo "#  gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-latest.tar.bz2.sig"; \
	 echo; \
	 echo "echo | gpg 2>/dev/null"; \
	 echo 'gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve $$0'; \
	 echo 'exit $$?' ; \
	 echo; \
	 gpg -ab -o - parallel-$(YYYYMMDD).tar.bz2; \
	) > parallel-latest.tar.bz2.sig

	gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-latest.tar.bz2.sig
	gpg --auto-key-locate keyserver --keyserver-options auto-key-retrieve parallel-$(YYYYMMDD).tar.bz2.sig

	../ftpsync/src/ftpsync parallel-$(YYYYMMDD).tar.bz2{,.sig,*asc} ftp://ftp-upload.gnu.org/incoming/alpha/

# This can take 7 minutes
	pushd /tmp; \
	rm -rf /tmp/parallel-$(YYYYMMDD)*; \
	while ! wget http://alpha.gnu.org/gnu/parallel/parallel-$(YYYYMMDD).tar.bz2 ; do sleep 2; done; \
	tar xjvf parallel-$(YYYYMMDD).tar.bz2; \
	cd parallel-$(YYYYMMDD); \
	./configure; \
	make -j && sudo make -j install; \
	pushd; \
	sudo cp /usr/local/bin/parallel /usr/local/bin/parallel-$(YYYYMMDD)

monitorman:
# If man page changed: open new pdfman
	inotifywait -qmre MOVED_TO -e CLOSE_WRITE --format %w%f . | parallel -uj1 'echo {=/\.pod$$/ or skip()=};make -j && sudo make install; pdfman {/.} &'

EXTRA_DIST = CITATION CREDITS
