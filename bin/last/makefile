CXXFLAGS = -msse4 -O3 -std=c++11 -pthread
all:
	@cd src && $(MAKE) CXXFLAGS="$(CXXFLAGS)"

prefix = /usr/local
exec_prefix = $(prefix)
bindir = $(exec_prefix)/bin
install: all
	mkdir -p $(bindir)
	cp bin/* $(bindir)

clean:
	@cd src && $(MAKE) clean

docs: doc/last-matrices.rst doc/last-seeds.rst

doc/last-matrices.rst: build/mat-doc.sh data/*.mat
	./build/mat-doc.sh data/*.mat > $@

doc/last-seeds.rst: build/seed-doc.sh data/*.seed
	cd data && ../build/seed-doc.sh [!R]*d *4-9* *8-10* *-11* *-12* > ../$@

tag:
	git tag -m "" `git rev-list HEAD | grep -c .`
