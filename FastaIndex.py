#/usr/bin/env python
"""FastA index (.fai) handler"""

class FastaIndex(object):
    """Facilitate Fasta index (.fai) operations compatible
    with samtools faidx (http://www.htslib.org/doc/faidx.html).
    """
    def __init__(handle):
        """ """
        self.handle = handle
        self.faidx  = handle.name + ".fai"
        if not os.path.isfile(self.faidx): # or older than fasta
            self._generate_index()
        else:
            self._load_fai()

    def _generate_index(self): 
        """Return fasta records"""
        header = ""
        seq = []
        self.id2stats = {}
        with open(self.faidx, 'w') as out:
            for l in iter(self.handle.readline, ''): 
                if l.startswith(">"):
                    if header:
                        self.id2stats[get_id(header)] = get_stats(header, seq, offset)
                        out.write("%s\t%s\n"%(get_id(header), "\t".join(map(str, stats))))
                    header = l
                    offset = handle.tell() 
                    seq = []
                else:
                    seq.append(l)

            if header: 
                self.id2stats[get_id(header)] = get_stats(header, seq, offset)
                out.write("%s\t%s\n"%(get_id(header), "\t".join(map(str, stats))))

    def _load_fai(self):
        """Load stats from faidx file"""
        self.id2stats = {}
        for l in open(self.faidx):
            ldata = l[:-1].split('\t')
            if len(ldata)<9:
                return {}
            rid = ldata[0]
            stats = map(int, ldata[1:])
            self.id2stats[rid] = stats

    def __getitem__(self, key):
        """x.__getitem__(y) <==> x[y]"""
        if key not in self.id2stats:
            raise KeyError
        # get offset info
        size, offset, linebases, linebytes = self.id2stats[key][:4]
        
        proxies = self._proxies
        if file_number in proxies:
            record = proxies[file_number].get(offset)
        else:
            if len(proxies) >= self._max_open:
                # Close an old handle...
                proxies.popitem()[1]._handle.close()
            # Open a new handle...
            proxy = self._proxy_factory(self._format, self._filenames[file_number])
            record = proxy.get(offset)
            proxies[file_number] = proxy
        if self._key_function:
            key2 = self._key_function(record.id)
        else:
            key2 = record.id
        if key != key2:
            raise ValueError("Key did not match (%s vs %s)" % (key, key2))
        return record


    def get_id(header):
        """Return seqid from header"""
        return header[1:].split()[0]

    def get_stats(header, seq, offset):
        """Compute stats"""
        errors = 0
        # get bases & bytes in line, ignoring last line
        if len(seq)>1:
            linebases = set(len(s.strip()) for s in seq[:-1])   
            linebytes = set(len(s) for s in seq[:-1])
            if len(linebases)>1:
                sys.stderr.write("[WARNING] Uneven line lengths in %s: %s\n"%(header, ",".join(map(str, linebases))))        
            linebases, linebytes = max(linebases), max(linebytes)
        else:
            linebases, linebytes = len(seq[0].strip()), len(seq[0])
        seq = "".join(s.strip() for s in seq)
        # count ACGT
        bases = {'A': 0, 'C': 0, 'G': 0, 'T': 0, 'N': 0}
        for b in seq.upper():
            if b in bases:
                try:
                    bases[b] += 1
                except:
                    errors += 1
        return (len(seq), offset, linebases, linebytes, \
                bases['A'], bases['C'], bases['G'], bases['T'])    
