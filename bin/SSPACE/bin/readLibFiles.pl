  #############################################################
  #Marten Boetzer 23-11-2011                                  #
  #SSPACE perl subscript readLibFiles.pl                      #
  #This script;                                               #
  #  -reads and converts original input sequences             #
  #############################################################

  use Storable;
  use File::Path;
  use File::Basename;
  use threads;

  my $seplines = ("-" x 60)."\n";

  my $libraryfile = $ARGV[0];
  my $base_name = $ARGV[1];
  my $thread = $ARGV[2];
  my $log = "$base_name/$base_name.logfile.txt";
  my $summaryfile = "$base_name/$base_name.summaryfile.txt";

  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";

#-------------------------------------------------LOOP THROUGH EACH LIBRARY IN LIBRARYFILE AND STORE AND FILTER READS
  open(FILELIB, "< $libraryfile");

  my ($library, $aligner, $fileA, $fileB, $insert_size, $insert_stdev, $reverse, $libResHash);
  my ($prevlibrary, $ctlib, $libcount) = ("",0, 0);
  &printMessage("\n=>".getDate().": Reading, filtering and converting input sequences of library file initiated\n");

  while(<FILELIB>){
    chomp;
    ($library, $aligner, $fileA, $fileB, $insert_size, $insert_stdev, $reverse) = split(/\s+/, $_);
    next if($library eq "");
    if($library ne $prevlibrary){
      $ctlib=0 if($prevlibrary ne "");
      push @libraries, $library;
    }
    $ctlib++;
    my $fname = "$base_name/reads/$base_name.$library.$aligner.file$ctlib";
    #Process multiple files at the same time if multithreaded option is set (-T parameter larger than 1)
    if($aligner ne "TAB"){
       my $thr = threads->create(\&readUnpairedFile, $fileA, $fname,$ctlib) if($library eq "unpaired");
       my $thr = threads->create(\&generateInputFiles, $library, $fileA, $fileB, $reverse, $fname, $ctlib, $aligner) if($library ne "unpaired");
       if(!(++$libcount % $thread)){
         foreach my $thr (threads->list()) {
           my @res = $thr->join();
           ($lib,$nreads) = split(/,/,$res[0]);
           $libResHash->{$lib}{'reads'}+=$nreads;
         }
       }
    }
    #if user has inserted a TAB file, calculate read statistics
    elsif($aligner eq "TAB"){
      open FILE, "$fileA" or die $!;
      my ($fileBaseName2, $dirName2, $fileExtension2) = fileparse($fileA);
      print "Reading tabfile: $fileBaseName2...\n";
      my $tabcount = 0;
      $tabcount++ while(<FILE>);
      $libResHash->{$library}{'reads'}+=$tabcount;
      close FILE;
    }
    $prevlibrary = $library;
  }
  #Process remaining threads
  foreach my $thr (threads->list()) {
      my @res = $thr->join();
      ($lib,$nreads) = split(/,/,$res[0]);
      $libResHash->{$lib}{'reads'}+=$nreads;
  }
  #Print read statistics to the summary file
  &printMessage("\n$seplines");
  foreach my $libs (@libraries){
    my $totcounter = $libResHash->{$libs}{'reads'};
    print SUMFILE "READING READS $libs:\n";
    print SUMFILE "$seplines\tTotal inserted reads = $totcounter \n$seplines\n" if($libs eq "unpaired");
    print SUMFILE "$seplines\tTotal inserted pairs = $totcounter \n$seplines\n" if($libs ne "unpaired");
  }
  close FILELIB;
  close SUMFILE;
  close LOG;

  mkpath('process_OK'); #make directory, indicating that process has run OK

#--------------------------------------------------

###CHECK FILE FORMAT: FASTA OR FASTQ? SEQUENCES ON ONE LINE OR MULTIPLE LINES?
sub generateInputFiles{
  my ($lib, $fileA, $fileB, $reverse, $fname, $libct, $aligner) = @_;
  my ($counterext, $Ncount, $countsinglet, $fastq, $step, $out) = (0,0,0,0,1000000, "");

  #check if file is zipped, otherwise unzip. C
  if ($fileA =~ /\.gz$/) {
    open(TEST, "gunzip -c  $fileA |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(TEST, "< $fileA");
  }
  #Check if file is fastQ or fastA
  $name = <TEST>; <TEST>; <TEST>; <TEST>;
  $fastq = 1 if ($name =~ /^[@]/);

  #check if sequences or on one line or multiline
  my ($ctlines, $cthead) = (0,0);
  while(<TEST>){
    my $head = $_;
    <TEST>;
    if($fastq){
      <TEST>;
      <TEST>;
    }
    $cthead++ if($head =~ /^@/ && $fastq);
    $cthead++ if($head =~ /^>/ && !$fastq);
    $ctlines++;
    last if($ctlines >= 1000);
  }
  close TEST;
  if($cthead != $ctlines){
    $out = getFastAMulti($fname, $fileA, $fileB, $lib, $base_name, $libct) if(!$fastq);
    $out = getFastQMulti($fname, $fileA, $fileB, $lib, $base_name, $libct) if($fastq);
  }else{
    $out = getFastXSingle($fname, $fileA, $fileB, $lib, $base_name, $libct, $fastq);
  }
  return $out;
}

###Sequence reads are on a single line, read and convert the pairs
sub getFastXSingle{
  my ($fname, $fileA, $fileB, $lib, $base_name, $libct, $fastq) = @_;
  my ($seq1, $seq2,$Ncount,$countsinglet, $step,$countSub) = ("","",0,0, 1000000,1);
  if ($fileA =~ /\.gz$/) {
    open(FILEA, "gunzip -c  $fileA |") || die "can't open .gz $fileA with gunzip";
  }else{
    open(FILEA, "< $fileA");
  }
  if ($fileB =~ /\.gz$/) {
    open(FILEB, "gunzip -c  $fileB |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(FILEB, "< $fileB");
  }
  open (OUTSINGLEFILE, ">$fname.1.fa") || die "Can't write to single file $fname.1.fa -- fatal\n";
  CounterPrint("Reading single line read-pairs $lib.$libct @ $countsinglet       ");
  while(<FILEA>) {
    <FILEB>;
    $seq1 = <FILEA>;
    $seq2 = <FILEB>;
    #FASTQ FORMAT
    <FILEA>,<FILEA>,<FILEB>,<FILEB> if ($fastq);
    print OUTSINGLEFILE ">read$countsinglet\n$seq1>read$countsinglet\n$seq2";
    if(++$countsinglet == $step){
      CounterPrint("Reading read-pairs $lib.$libct @ $countsinglet         ");
      $step = $step + 1000000;
      close OUTSINGLEFILE;
      $countSub++;
      open (OUTSINGLEFILE, ">$fname.$countSub.fa") || die "Can't write to single file $fname.$countSub.fa -- fatal\n";
    }

  }
  CounterPrint("\n") if($thread <= 1);
  CounterPrint((" " x 40));
  close OUTSINGLEFILE;
  close FILEB;
  close FILEA;
  return "$lib,$countsinglet";
}

###read-in and convert fastA formatted paired-reads where sequences are on multiple lines
sub getFastAMulti{
  my ($fname,$fileA, $fileB, $lib, $base_name, $libct) = @_;
  my ($countSub, $seq1, $seq2, $line1, $Ncount,$countsinglet, $step, $res1, $res2) = (1, "","",0,0, 1000000);
  if ($fileA =~ /\.gz$/) {
    open(FILEA, "gunzip -c  $fileA |") || die "can't open .gz $fileA with gunzip";
  }else{
    open(FILEA, "< $fileA");
  }
  if ($fileB =~ /\.gz$/) {
    open(FILEB, "gunzip -c  $fileB |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(FILEB, "< $fileB");
  }
  open (OUTSINGLEFILE, ">$fname.1.fa") || die "Can't write to single file $fname.1.fa -- fatal\n";
  CounterPrint("Reading multiline read-pairs $lib.$libct @ $countsinglet         ");
  while($line1 = <FILEA>) {
    chomp $line1;
    if($line1 =~ /^>/ || eof(FILEA)){
      $seq1.=$line1 if(eof(FILEA));
      if($seq1 ne ""){
        my ($line2, $seq2) = ("","");
        while($line2 = <FILEB>) {
            chomp $line2;
            if($line2 =~ /^>/ || eof(FILEB)){
              $seq2.=$line2 if(eof(FILEB));
              if($seq2 ne ""){
                if(++$countsinglet == $step){
                  CounterPrint("Reading multiline read-pairs $lib.$libct @ $countsinglet         ");
                  $step = $step + 1000000;
                  close OUTSINGLEFILE;
                  $countSub++;
                  open (OUTSINGLEFILE, ">$fname.$countSub.fa") || die "Can't write to single file $fname.$countSub.fa -- fatal\n";
                }
                print OUTSINGLEFILE ">read$countsinglet\n$seq1\n>read$countsinglet\n$seq2\n";
                last;
              }
            $seq2 = "";
             }
          else{
            $seq2.=$line2;
          }
        }
      }
      $seq1 = "";
    }
    else{
      $seq1.=$line1;
    }
  }
  close FILEB;
  close FILEA;
  return "$lib,$countsinglet";
}

###read-in and convert fastQ formatted paired-reads where sequences are on multiple lines
sub getFastQMulti{
  my ($fname, $fileA, $fileB, $lib, $base_name, $libct) = @_;
  my ($countSub, $seq1, $seq2, $line1, $Ncount,$countsinglet, $step, $res1, $res2) = (1, "","","",0,0, 1000000);
  if ($fileA =~ /\.gz$/) {
    open(FILEA, "gunzip -c  $fileA |") || die "can't open .gz $fileA with gunzip";
  }else{
    open(FILEA, "< $fileA");
  }
  if ($fileB =~ /\.gz$/) {
    open(FILEB, "gunzip -c  $fileB |") || die "can't open .gz $fileB with gunzip";
  }else{
    open(FILEB, "< $fileB");
  }
  open (OUTSINGLEFILE, ">$fname.1.fa") || die "Can't write to single file $fname.1.fa -- fatal\n";
  CounterPrint("Reading multiline read-pairs $lib.$libct @ $countsinglet         ");
  while($line1 = <FILEA>) {
    chomp $line1;
    if($line1 =~ /^\+/){
      if($seq1 ne ""){
        my ($line2, $seq2) = ("","");
        while($line2 = <FILEB>) {
            chomp $line2;
            if($line2 =~ /^\+/){
              if($seq2 ne ""){
                if(++$countsinglet == $step){
                  CounterPrint("Reading multiline read-pairs $lib.$libct @ $countsinglet         ");
                  $step = $step + 1000000;
                  close OUTSINGLEFILE;
                  $countSub++;
                  open (OUTSINGLEFILE, ">$fname.$countSub.fa") || die "Can't write to single file $fname.$countSub.fa -- fatal\n";
                }
                print OUTSINGLEFILE ">read$countsinglet\n$seq1\n>read$countsinglet\n$seq2\n";
                last;
              }
            $seq2 = "";
          }elsif($line2 =~ /^@/){
            $seq2 = "";
          }
          else{
            $seq2.=$line2;
          }
        }
      }
      $seq1 = "";
    }elsif($line1 =~ /^@/){
      $seq1 = "";
    }
    else{
      $seq1.=$line1;
    }
  }
  close FILEB;
  close FILEA;
  return "$lib,$countsinglet";
}



#------------------READ UNPAIRED SINGLE READS FILES FOR EXTENSION

sub readUnpairedFile{
  my ($file, $fname, $libct) = @_;
  if ($file =~ /\.gz$/) {
    open(INUNPAIRED, "gunzip -c  $file |") || die "can't open .gz $file with gunzip";
  }else{
    open(INUNPAIRED, "< $file") || die "Can't open $file -- fatal\n";
  }
  open OUTFILEExt, "> $fname.1.fa";

  my ($seq1, $counter, $step, $fastq) = ("",0, 2000000,0);

  open(TEST, "< $file");
  my $name = <TEST>;
  close TEST;
  $fastq = 1 if ($name =~ /^[@]/);
  while(<INUNPAIRED>) {
    $seq1 = <INUNPAIRED>;

    #FASTQ FORMAT
    if ($fastq){
      <INUNPAIRED>; <INUNPAIRED>;
    }
    # ELSE FASTA FORMAT
    print OUTFILEExt ">$counter\n$seq1";
    if(++$counter == $step){
      CounterPrint("Reading unpaired reads @ $counter                     ");
      $step = $step + 2000000;
      close OUTFILEExt;
      $countSub++;
      open (OUTFILEExt, ">$fname.$countSub.fa") || die "Can't write to single file $fname.$countSub.fa -- fatal\n";
    }
  }
  CounterPrint("                ");
  close OUTFILEext;
  close INUNPAIRED;
  
  return "unpaired,$counter";
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

#########END readLibFiles.pl