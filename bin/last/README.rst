LAST: find & align related regions of sequences
===============================================

LAST is designed for moderately large data (e.g. genomes, DNA reads,
proteomes).  It's especially good at:

* Finding rearrangements and recombinations: we believe last-split_
  does that more rigorously than anything else.

* Finding DNA-versus-protein related regions, especially protein_
  fossils_.

* Unusual data, e.g. AT-rich DNA, because we can fit_ parameters to
  the data and calculate significance_.

* Sensitive DNA-DNA search, due to fitting_, sensitive seeding_, and
  calculating significance_.

It can also: indicate the confidence/uncertainty of each column in an
alignment, and use sequence quality data in a rigorous fashion.

Usage
-----

Please see the cookbook_.  **Warning:** this documentation may not
apply to older versions of LAST!  You can see your version with::

  lastal --version

Install
-------

You can install it from bioconda_ or `Debian Med`_, or like this...

Download the highest version number from
https://gitlab.com/mcfrith/last/-/tags.  Using the command line, go
into the downloaded directory and type::

  make

This assumes you have a C++ compiler.  On Linux, you might need to
install a package called "g++".  On Mac, you might need to install
command-line developer tools.  On Windows, you might need to install
Cygwin.  You might also need to install something like "zlib-devel".

For ARM CPUs, the default "make" seems to work in some cases but not
others (sigh).  This seems to be good for ARM::

  make CXXFLAGS="-mcpu=native -O3 -pthread"

It's possible to specify a compiler like this: ``make CXX=MyOtherCompiler``.
If you re-run ``make`` in different ways, it may be good to do ``make clean``
first, to remove any previously-made files.

The programs are in the ``bin`` directory.  For convenient usage, set
up your computer to find them automatically.  Some possible ways:

* Copy the programs to a standard directory: ``sudo make install``
  (using "sudo" to request administrator permissions).

* Copy the programs to your personal bin directory: ``make install prefix=~``

* Adjust your `PATH variable`_.

You might have to log out and back in before your computer recognizes
the new programs.

Further info
------------

Details & citation: `LAST papers`_

LAST is distributed under the GNU General Public License, either
version 3 of the License, or (at your option) any later version.

LAST is brought to you by:

* `Computational Omics Research Team`_, AIRC_
* GSFS_, `University of Tokyo`_
* `AIST-Waseda University CBBD-OIL`_

.. _fit:
.. _fitting: doc/last-train.rst
.. _last-split: doc/last-split.rst
.. _seeding: doc/last-seeds.rst
.. _significance: doc/last-evalues.rst
.. _cookbook: doc/last-cookbook.rst
.. _LAST papers: doc/last-papers.rst
.. _protein: https://doi.org/10.1109/TCBB.2022.3177855
.. _fossils: https://doi.org/10.1093/molbev/msac068
.. _bioconda: https://bioconda.github.io/
.. _Debian Med: https://www.debian.org/devel/debian-med/
.. _PATH variable: https://en.wikipedia.org/wiki/PATH_(variable)
.. _Computational Omics Research Team: https://www.airc.aist.go.jp/en/cort/
.. _AIRC: https://www.airc.aist.go.jp/en/
.. _GSFS: https://www.k.u-tokyo.ac.jp/index.html.en
.. _University of Tokyo: https://www.u-tokyo.ac.jp/en/
.. _AIST-Waseda University CBBD-OIL: https://unit.aist.go.jp/cbbd-oil/en/
