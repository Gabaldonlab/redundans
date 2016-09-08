
  ###############################################################
  #Marten Boetzer 23-11-2010                                    #
  #SSPACE perl subscript ExtendOrFormatContigs.pl               #
  #This script, based on the the -x parameter;                  #
  #  -Formats the contigs to appropriate format (-x 0)          #
  #  -Extends the contigs with available unmapped reads (-x 1)  #
  #If -z option is set, contigs below -z are filtered out       #
  #                                                             #
  ###############################################################

  use strict;
  use File::Basename;
  use File::Path;
  use threads;
  use threads::shared;
  my $countReads :shared;

  my $seplines = ("-" x 60)."\n";

#get input parameters
  my $contig = $ARGV[0];
  my $base_name = $ARGV[1];
  my $extending = $ARGV[2];
  my $filecontig = $ARGV[3];
  my $min_coverage = $ARGV[4];
  my $min_overlap = $ARGV[5];
  my $min_base_ratio = $ARGV[6];
  my $Bin = $ARGV[7];
  my $minContigLength = $ARGV[8];
  my $gaps = $ARGV[9];
  my $threads = $ARGV[10];
  my ($bin);

#update log and summary files
  my $log = "$base_name/$base_name.logfile.txt";
  my $summaryfile = "$base_name/$base_name.summaryfile.txt";

  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open (LOG, ">>$log") || die "Can't write to logfile$log -- fatal\n";
  my $filenameOutExt = $base_name . ".singlereads.fasta";

#extend contigs or format contigs
  if($extending == 1){
    &ExtendContigs($base_name, $filecontig, $filenameOutExt, $Bin);
    print SUMFILE "\n" if($minContigLength > 0);
    &FormatContigs() if($minContigLength > 0);
  }else{
    &FormatContigs();
  }

  close SUMFILE;
  close LOG;
  
  mkpath('process_OK');
#--------------------------------------------------

###EXTEND CONTIGS WITH UNMAPPED READS
sub ExtendContigs{
  my ($base_name, $filecontig, $filenameOutExt, $Bin) = @_;
  #-------------------------------------------------NOW MAP SINGLE READS TO INITIAL CONTIGS FILE.
  &printMessage("\n=>".getDate().": Extending of contigs with unmapped reads\n");
  my $readfile = "$base_name/reads/" . $filenameOutExt;
  &getUnmappedReads($filecontig, $readfile);
  #-------------------------------------------------CONTIG EXTENSION USING UNMAPPED PAIRS STORED IN $SET
  &printMessage("\n=>".getDate().": Contig extension initiated\n");
  my $outfileTig =  "$base_name/intermediate_results/" . $base_name .  ".extendedcontigs.fasta";
  my $outfileEvidence =  "$base_name/intermediate_results/" . $base_name .  ".extension_evidence.txt";

  open (TIG, ">$outfileTig") || die "Can't write to $outfileTig -- fatal\n";
  open (EVI, ">$outfileEvidence") || die "Can't write to $outfileEvidence -- fatal\n";
  #--------------------------------------------ASSEMBLY START

  #obtain the contig sequences and extend them unmapped reads
  ASSEMBLY:
    open(IN,$filecontig) || die "Can't open $filecontig -- fatal\n";
    my ($seq,$exttig_count, $counter, $extend_bases, $orig_mer, $prevhead, $step) = ("",0,0,0,0,'',100);
    while(<IN>){
      s/\r\n/\n/;
      chomp;
      $seq.= uc($_) if(eof(IN));
      if (/\>(\S+)/ || eof(IN)){
         my $head=$1;
         $orig_mer = length($seq);
         if($seq ne ''){
           print EVI ">contig$counter|orig:$prevhead\n";
           ($seq) = doExtension("3", $orig_mer, $seq, $min_overlap, $min_coverage, $min_base_ratio, $counter) if($orig_mer >= $min_overlap);
           
           my $seqrc = reverseComplement($seq);
           ($seqrc) = doExtension("5", $orig_mer, $seqrc, $min_overlap, $min_coverage, $min_base_ratio, $counter) if($orig_mer >= $min_overlap);
           
           my $leng = length($seqrc);
           my $reversetig = reverseComplement($seqrc);                   ### return orientation of the sequence, as inputted
           if($leng > $orig_mer){ #print contigs to file if contig is extended
             $extend_bases+=($leng-$orig_mer);
             printf TIG ">extcontig%i|size%i|prevsize%i|seed:$prevhead\n%s\n", ($counter,$leng,$orig_mer,$reversetig);
             $exttig_count++;
           }else{
             printf TIG ">contig%i|size%i|seed:$prevhead\n%s\n", ($counter,$leng,$reversetig);    #print singlets to file
           }
         }
         if(++$counter == $step){
            CounterPrint($counter);
            $step = $step + 100;
         }
         $prevhead = $head;
         $seq='';
      }else{
         $seq .= uc($_);
      }
   }
  CounterPrint("                ");
  print SUMFILE "CONTIG EXTENSION:\n$seplines\n";
  print SUMFILE "\tNumber of contig sequences =".($counter-1). "\n";
  print SUMFILE "\tNumber of unmapped reads used for contig extension = $countReads\n";
  print SUMFILE "\tNumber of contigs extended = $exttig_count\n";
  print SUMFILE "\tNumber of bases extended = $extend_bases\n".$seplines;
  close IN;
  close EVI;
  $filecontig = $outfileTig;
  if($@){
     my $message = $@;
     &printMessage("\nSomething went wrong running $0 ".getDate()."\n$message\n");
  }
  close EVI;
  close TIG;
}

###STORE CONTIGS TO APPROPRIATE FORMAT WHEN CONTIGS WILL NOT BE EXTENDED
sub FormatContigs{
   &printMessage("\n=>".getDate().": Storing contigs to format for scaffolding\n");
   open (TIG, ">$contig") || die "Can't write to $contig -- fatal\n";
   die "Can't read file $filecontig. Check your permissions\n" if(!-r $filecontig);
   open(IN,"$filecontig") || die "Can't open $filecontig -- fatal\n";
   my ($counter, $seq, $prevhead, $step) = (0,'','', 100);
   while(<IN>){
      s/\r//;
      chomp;
      $seq.= uc($_) if(eof(IN));
      if (/\>(\S+)/ || eof(IN)){
        my $head=$1;
        my $length_seq = length($seq);
        if($seq ne '' && $length_seq >= $minContigLength){    #if length of contig meets the minimum contig length (-z option)
          if(++$counter == $step){
            CounterPrint($counter);
            $step = $step + 100;
          }
          printf TIG ">contig%i|size%i|seed:$prevhead\n%s\n", ($counter,$length_seq,$seq);
        }
        $prevhead = $head;
        $seq = '';
      }else{
         $seq .= uc($_);
      }
   }
   CounterPrint("                ");
   close IN;
   close TIG;
}

###EXTEND CONTIGS
sub doExtension{
  my ($direction, $orig_mer, $seq, $min_overlap, $min_coverage, $min_base_ratio, $tig_count) = @_;
  my ($previous, $extended) = ($seq,1);
  
  CONSENSUS:
  while($extended){
    $extended=0;
    my $subseq = substr($seq, -$min_overlap);
    my $revseq = reverseComplement($subseq);
    my $overhang;
    #get one basepair overhang of the contig and determine how many times each nucleotide occurs
    $overhang->{'A'} = $bin->{$subseq."A"}+$bin->{"T$revseq"};
    $overhang->{'C'} = $bin->{$subseq."C"}+$bin->{"G$revseq"};
    $overhang->{'G'} = $bin->{$subseq."G"}+$bin->{"C$revseq"};
    $overhang->{'T'} = $bin->{$subseq."T"}+$bin->{"A$revseq"};
    
    my $coverage = $overhang->{'A'}+$overhang->{'C'}+$overhang->{'G'}+$overhang->{'T'};
    print EVI "\tdir$direction: Total:$coverage\tA:$overhang->{'A'}\tT:$overhang->{'T'}\tG:$overhang->{'G'}\tC:$overhang->{'C'}";
    my ($ct_dna, $previous_bz, $extend_nuc) = (0, "", "");
    #get the two best nucleotides, and see if the best nucleotide meets the criteria (above minimal coverage and below minimum base ratio between the two best bases
    BASE:
    foreach my $bz (sort {$overhang->{$b}<=>$overhang->{$a}} keys %$overhang){
      if($ct_dna == 1){## the two most abundant bases at that position
        my $bestcoverage = $overhang->{$previous_bz} + $overhang->{$bz};
        if($overhang->{$previous_bz} < $min_coverage){      #does the best nucleotide occur more than the minimum coverage
          $extended = 0;
          print EVI " => low coverage\n";
          last CONSENSUS;
        }
        if(($overhang->{$previous_bz} / $bestcoverage) >= $min_base_ratio){#determine the ratio between the two best bases
          $extend_nuc = "$previous_bz";
          print EVI " => extending with $extend_nuc\n";
          last BASE;
        }
        print EVI " => below ratio\n";
        last CONSENSUS;
      }
      $previous_bz = $bz;
      $ct_dna++;
    }
    deleteData($subseq."$extend_nuc");
    $seq = $seq . $extend_nuc;
    $extended = 1;
  }###while get the OK for extension
  return $seq;
}


###DELETE READ DATA IF IT HAS BEEN USED FOR EXTENDING A CONTIG
sub deleteData {
   my ($sequence) = @_;
   my $comp_seq = reverseComplement($sequence);
   #remove k-mer and its reverse complement from hash table
   delete $bin->{$comp_seq};
   delete $bin->{$sequence};
}

#obtain all bowtie and bwa readfiles, map them to the contigs and get unmapped reads
sub getUnmappedReads{
  my ($contigFile, $readfiles) = @_;
  #determine number of readfiles that should be mapped with bowtie

  my @filesbowtie = <$base_name/reads/$base_name.*.bowtie*fa>;
  my $ctlib = 0;
  if($#filesbowtie >= 0){
    my $bowtieout = $base_name . ".bowtieIndexExt";
    my $bowtiepath = "$Bin"."/bowtie/bowtie";
    $bowtiepath =~ s/ /\\ /g;
    my $bowbuildpath = $bowtiepath."-build";
    die "Contig file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));
    &printMessage("\n=>".getDate().": Building Bowtie index for contigs\n");
     system("$bowbuildpath $contigFile $base_name/alignoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?"; # returns exit status values
    &printMessage("\n=>".getDate().": Mapping reads to contigs with Bowtie\n");
    #for each readfile, map the reads to the contigs and obtain unmapped reads with function readSam
    foreach my $fname (@filesbowtie){
      my $procline = "$bowtiepath -p 1 -v $gaps -m 1 $base_name/alignoutput/$bowtieout -f $fname --quiet --refidx -S --sam-nohead|";
      my $thr = threads->create(\&readSam, "", $procline);
      if(!(++$ctlib % $threads)){
        readThreadOutput();
      }
    }
  }
  #determine number of reads that should be mapped with bwa (either bwasw or normal bwa)
  my @filesbwa = <$base_name/reads/$base_name.*.bwa.*fa>;
  my @filesbwasw = <$base_name/reads/$base_name.*.bwasw.*fa>;
  if($#filesbwa >= 0 || $#filesbwasw >= 0){
    my $bwapath = "$Bin"."/bwa/bwa";
    $bwapath =~ s/ /\\ /g;
    my $bwaout = $base_name . ".bwaIndexExt";
    die "Contig file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));
    &printMessage("\n=>".getDate().": Building BWA index for contigs\n");
    my $filesize = -s "$contigFile";
    my $index = "bwtsw";
    $index = "is" if($filesize <= 10000000);
    mkpath("$base_name/tmp.$base_name");
    open(STDERR, ">$base_name/tmp.$base_name/tmpbwa_logfile_extension");
    system("$bwapath index -a $index $contigFile -p $base_name/alignoutput/$bwaout") == 0 || die "\nBwa error; $?"; # returns exit status values
    &printMessage("\n=>".getDate().": Mapping reads to contigs with BWA\n");
    #for each readfile, map the reads to the contigs and obtain unmapped reads with function readSam
    foreach my $fname (@filesbwa){
      my $bwaoutputaln = "$base_name/alignoutput/$base_name.extension.$ctlib.bwa";
      my $procline = "$bwapath aln $base_name/alignoutput/$bwaout $fname > $bwaoutputaln";# >/dev/null 2>&1";
      my $samseline = "$bwapath samse $base_name/alignoutput/$bwaout $bwaoutputaln $fname |";
      my $thr = threads->create(\&readSam, $procline, $samseline);
      if(!(++$ctlib % $threads)){
        readThreadOutput();
      }
    }
    #for each readfile, map the reads to the contigs and obtain unmapped reads with function readSam
    foreach my $fname (@filesbwasw){
      my $procline = "$bwapath bwasw $base_name/alignoutput/$bwaout $fname |";
      my $thr = threads->create(\&readSam, "", $procline);
      if(!(++$ctlib % $threads)){
        readThreadOutput();
      }
    }
  }
  close STDERR;
  close FILELIB;
  readThreadOutput();
}

#read sam alignment output and find reads that do not match to any contig
sub readSam{
  my ($aligninput, $samse) = @_;
  my ($subbin, $tig_length);
  system("$aligninput") == 0 || die "\nBwa error; $?" if($aligninput ne "");
  open(SAM, "$samse") || die "Can't open bwa output -- fatal\n";
  #go through mapping results
  while(my $line = <SAM>){
    my @t = split(/\t/, $line);
    next if($_ =~ /^@/);
    if($t[2] eq "*"){
      $subbin->{$t[9]}++;
      $countReads++;
    }
  }
  close SAM;
  return $subbin;
}

#read the found unmapped reads of each thread and store it in memory
sub readThreadOutput{
  foreach my $thr (threads->list()) {
    my @res = $thr->join();
    my $partHash = $res[0];
    foreach my $seq (keys %$partHash){
      my $subct=0;
      while($subct < length($seq)-$min_overlap){
        my $subnor = substr($seq, $subct, $min_overlap+1);
        $subct++;
        next if(index($subnor,"N") >= 0);
        if(defined $bin->{$subnor}){
          $bin->{$subnor}++;
        }else{
          my $subrv =reverseComplement($subnor);
          $bin->{$subrv}++;
        }
      }
    }
  }
}

###FUNCTION TO REVERSE COMPLEMENT A SEQUENCE
sub reverseComplement{
   $_ = shift;
   tr/ATGC/TACG/;
   return (reverse());
}

###PRINTS A COUNTER ON THE SCREEN AND OVERWRITES PREVIOUS LINE
sub CounterPrint{
  my $countingMessager = shift;
  print "\r$countingMessager";
  $|++;
}

###FUNCTION TO PRINT MESSAGES TO THE SCREEN AND TO THE LOG FILE
sub printMessage{
  my $message = shift;
  print $message;
  print LOG $message;
}

###FUNCTION TO GET THE CURRENT DATE
sub getDate{
  my $date = scalar(localtime);
  return $date;
}

###FLUSHES THE SUMMARY AND LOG FILE
sub FlushFiles{
  select((select(SUMFILE), $| = 1)[0]);
  select((select(LOG), $| = 1)[0]);
  $|++;
}

sub checkStatus{
  &printMessage(("*" x 50)."\n\nProcess failed on ".getDate()."\n\n\n"), exit if(!(-d "process_OK"));
  rmtree(["process_OK", 'blurfl/quux']);
}

#########END ExtendOrFormatContigs.pl