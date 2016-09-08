  ###################################################
  #Marten Boetzer 6-11-2012                         #
  #SSPACE perl subscript PairingAndScaffolding.pl   #
  #This script;                                     #
  #  -reads the contig sequences in a hash          #
  #  -stores Bowtie/BWA output in a hash            #
  #  -pairs the contigs                             #
  #  -generates scaffolds                           #
  ###################################################

  use strict;
  use Storable;
  use File::Path;
  use File::Basename;
  use Cwd;
  use threads;
  use threads::shared;
  my $totalReads :shared;
  my $totalPairsFound :shared;

  my $cwd = getcwd;
  my $seplines = ("-" x 60)."\n";

  my $Bin = $ARGV[0];
  my $gaps = $ARGV[1];
  my $contig = $ARGV[2];
  my $base_name = $ARGV[3];
  my $verbose = $ARGV[4];
  my $library = $ARGV[5];
  my $insert_size = $ARGV[6];
  my $min_allowed = $ARGV[7];
  my $scaffold = $ARGV[8];
  my $min_links = $ARGV[9];
  my $max_link_ratio = $ARGV[10];
  my $ori = $ARGV[11];
  my $threads = $ARGV[12];
  my $tab = $ARGV[13];
  my $tabfile = $ARGV[14];
  my $origctg = $ARGV[15];
  my $prev_evidence = $ARGV[16];

  my ($low_iz, $up_iz) = ($insert_size + $min_allowed, $insert_size - $min_allowed);
  my $subContigLim = ($up_iz * 2) + 3;
  my $lenlimit = (($up_iz * 2));
  my $bowtiefile = "$base_name/alignoutput/" . $base_name . ".$library.mapped";
  my $log = "$base_name/$base_name.logfile.txt";
  my $summaryfile = "$base_name/$base_name.summaryfile.txt";
  my $issues = "$base_name/pairinfo/" . $base_name . ".$library.pairing_issues";
  my $distribution = "$base_name/pairinfo/" . $base_name . ".$library.pairing_distribution.csv";

  my ($total_for_median, $ct_illogical, $ct_ok_contig, $ct_ok_pairs, $ct_problem_pairs, $ct_iz_issues)= (0,0,0,0,0,0);
  my ($pair,$err,$track_insert, $tigOnScafHash, $tigHash, $count);
  my ($totalPairsUsed, $pair_found) = (0,0,0,0);

  open (LOG, ">>$log") || die "Can't write to $log -- fatal\n";
  open (SUMFILE, ">>$summaryfile") || die "Can't open $summaryfile -- fatal\n";
  open(PET, ">$issues") || die "Can't open $issues for writing -- fatal\n";
#-------------------------------------------------READ CONTIGS INTO HASH AND STORE THEIR LENGTH. NEXT; PAIR THE CONTIGS

  #if a tab file is inserted with pre-defined mapping results, update the contigs if the original contigs are within scaffolds
  if($tab){
    parseEvidenceFile($prev_evidence) if($prev_evidence ne '');
    &updateContigs($origctg);
  }
  #determine lengths of the contigs/scaffolds
  my ($contigstored, $tig_length) = &readFileContigHash($contig);
  #find all bowtie/bwa files of the current library and map them to the contigs/scaffolds
  if(!$tab){
    my $contigFile = processContig($contig) ;
    my ($ctlib,$libnum)=(0,1);
    #find all bowtie files
    my @filesbowtie = <$base_name/reads/$base_name.$library.bowtie*fa>;
    if($#filesbowtie >= 0){
      my $bowtieout = $base_name . ".$library.bowtieIndex";
      my $bowtiepath = "$Bin"."/bowtie/bowtie";
      $bowtiepath =~ s/ /\\ /g;
      my $bowbuildpath = $bowtiepath."-build";
      die "Contig file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));
      &printMessage("\n=>".getDate().": Building Bowtie index for contigs\n");
      system("$bowbuildpath $contigFile $base_name/alignoutput/$bowtieout --quiet --noref") == 0 || die "\nBowtie-build error; $?"; # returns exit status values
      &printMessage("\n=>".getDate().": Mapping reads to contigs with Bowtie\n");
      foreach my $fname (@filesbowtie){
        my $procline = "$bowtiepath -p 1 -v $gaps -m 1 $base_name/alignoutput/$bowtieout --suppress 6,7 -f $fname --quiet --refidx |";
        my $thr = threads->create(\&readBowtie, $procline);
        if(!(++$ctlib % $threads)){
          readThreadOutput();
        }
      }
    }
    #find all bwa/bwasw files
    my @filesbwa = <$base_name/reads/$base_name.$library.bwa.*fa>;
    my @filesbwasw = <$base_name/reads/$base_name.$library.bwasw.*fa>;
    if($#filesbwa >= 0 || $#filesbwasw >= 0){
      my $bwapath = "$Bin"."/bwa/bwa";
      $bwapath =~ s/ /\\ /g;
      my $bwaout = $base_name . ".$library.bwaIndex";
      die "Contig file ($contigFile) not found. Exiting...\n" if(!(-e $contigFile));
      &printMessage("\n=>".getDate().": Building BWA index for contigs\n");
      my $filesize = -s "$contigFile";
      my $index = "bwtsw";
      $index = "is" if($filesize <= 10000000);
      open(STDERR, ">$base_name/tmp.$base_name/tmpbwa_logfile");
      system("$bwapath index -a $index $contigFile -p $base_name/alignoutput/$bwaout") == 0 || die "\nBwa error; $?"; # returns exit status values
      &printMessage("\n=>".getDate().": Mapping reads to contigs with BWA\n");
      foreach my $fname (@filesbwa){
        my $bwaoutputaln = "$base_name/alignoutput/$base_name.$library.$ctlib.bwa";
        my $procline = "$bwapath aln $base_name/alignoutput/$bwaout $fname > $bwaoutputaln";# >/dev/null 2>&1";
        my $samseline = "$bwapath samse $base_name/alignoutput/$bwaout $bwaoutputaln $fname |";
        my $thr = threads->create(\&readBwa, $procline, $samseline);
        if(!(++$ctlib % $threads)){
          readThreadOutput();
        }
      }
      foreach my $fname (@filesbwasw){
        my $procline = "$bwapath bwasw $base_name/alignoutput/$bwaout $fname |";
        my $thr = threads->create(\&readBwa, "", $procline);
        if(!(++$ctlib % $threads)){
          readThreadOutput();
        }
      }
    }
    readThreadOutput();
    print SUMFILE "\nMAPPING READS TO CONTIGS:\n";
    print SUMFILE "$seplines\tNumber of single reads found on contigs = $totalReads\n";
    my $read_number_message = "\tNumber of read-pairs used for pairing contigs / total pairs = $totalPairsUsed / $totalPairsFound\n";
    printf SUMFILE $read_number_message.$seplines."\n";
  }else{
    parseTabFile($tigHash);
  }
  &printResultsPairing();

#-------------------------------------------------BUILDING SCAFFOLDS
  &buildScaffolds($pair, $tig_length, $verbose, $scaffold, $library);
  ($pair, $tig_length) = ('',''); undef $pair; undef $tig_length;

  close SUMFILE;
  close LOG;   
  mkpath('process_OK');

#-------------------------------------------------
###
sub readThreadOutput{
  foreach my $thr (threads->list()) {
    my @res = $thr->join();
    my $partHash = $res[0];
     foreach my $h (keys %$partHash){
       my ($track1, $seq1, $track2, $seq2) = split(",", $h);
     #  my $combined = "$seq1:$seq2";
     #  my $revcombined = reverseComplement($combined);
       #if(!$count->{$combined} && !$count->{$revcombined}){
      #   $count->{$combined}++;
         pairContigs($track1, $track2, $seq1, $seq2);
         $totalPairsUsed++;
      # }
     }
  }
}

###FUNCTION TO PARSE THE TAB FILE
sub parseTabFile{
  my ($tigHash) = @_;
  open(TAB, "$tabfile") || die "Can't open $tabfile for reading -- fatal\n";
  my ($ct_both,$step) = (0,1000000);
  &printMessage("\n=>".getDate().": Parsing Tab file\n");
  while(<TAB>){
    chomp;
    if(++$ct_both == $step){
      CounterPrint($ct_both);
      $step = $step + 1000000;
    }
    my ($tig1,$start1,$end1,$tig2,$start2,$end2) = split(/\t/);

    #check if the contig in the tab file is also present in the inserted contig fasta file
    if(!defined($tigHash->{$tig1})){
      die "\nERROR: could not find an header containing $tig1 at line number $ct_both. Exit...\n";
    }
    if(!defined($tigHash->{$tig2})){
      die "ERROR: could not find an header containing $tig2 at line number $ct_both. Exit...\n";
    }
    my $ctg1 = $tigHash->{$tig1};
    my $ctg2 = $tigHash->{$tig2};
    my ($track1, $track2) =  ("","");
    
    #if multiple libraries were used, update the contig positions in the TAB File by finding its position in the scaffolds
    if($prev_evidence ne ''){
      $start1 = $start1 + $tigOnScafHash->{$ctg1}{'begin'};
      $end1 = $end1 + $tigOnScafHash->{$ctg1}{'begin'};
      $start2 = $start2 + $tigOnScafHash->{$ctg2}{'begin'};
      $end2 = $end2 + $tigOnScafHash->{$ctg2}{'begin'};
      if($tigOnScafHash->{$ctg1}{'direction'} eq "r"){
        my $tmp_start = ($tigOnScafHash->{$ctg1}{'end'}+$tigOnScafHash->{$ctg1}{'begin'}) - $start1;
        my $tmp_end = ($tigOnScafHash->{$ctg1}{'end'}+$tigOnScafHash->{$ctg1}{'begin'}) - $end1;
        $start1 = $tmp_start;
        $end1 = $tmp_end;
      }
      if($tigOnScafHash->{$ctg2}{'direction'} eq "r"){
        my $tmp_start = ($tigOnScafHash->{$ctg2}{'end'}+$tigOnScafHash->{$ctg2}{'begin'}) - $start2;
        my $tmp_end = ($tigOnScafHash->{$ctg2}{'end'}+$tigOnScafHash->{$ctg2}{'begin'}) - $end2;
        $start2 = $tmp_start;
        $end2 = $tmp_end;
      }

      $ctg1 = $tigOnScafHash->{$ctg1}{'scaf'};
      $ctg2 = $tigOnScafHash->{$ctg2}{'scaf'};
     if($start1 < $up_iz || ($end1 > ($tig_length->{$ctg1}-$up_iz))){
       $track1 = "$ctg1"."|$start1"."|$end1";
     }
     if($start2 < $up_iz || ($end2 > ($tig_length->{$ctg2}-$up_iz))){
       $track2 = "$ctg2"."|$start2"."|$end2";
     }
    }else{ #if it is the first library, just use the positions in the TAB file
      if($start1 < $up_iz || ($end1 > ($tig_length->{$ctg1}-$up_iz))){
        $track1 = "$ctg1"."|$start1"."|$end1";
      }
      if($start2 < $up_iz || ($end2 > ($tig_length->{$ctg2}-$up_iz))){
        $track2 = "$ctg2"."|$start2"."|$end2";
      }

    }
    #pair the contigs based on the information provided in the TAB file
    pairContigs($track1, $track2, "seq$ct_both.1", "seq$ct_both.2") if($track1 ne "" && $track2 ne "");
  }
  CounterPrint("                ");
}

###FUNCTION TO STORE ONLY THE EDGES OF THE CONTIGS. ONLY READS ARE MAPPED TO THESE EDGES, SAVING TIME FOR BUILDING THE INDEX WITH BOWTIE, AND MAPPING THE READS TO THE CONTIGS
sub processContig{
  my ($contigfile) = @_;
  open(IN,$contigfile) || die "can't read $contigfile -- fatal\n";
  my $contigfilesub = "$base_name/tmp.$base_name/subset_contigs.fasta";
  open(OUT,">$contigfilesub") || die "can't write to $contigfilesub -- fatal\n";
  my ($seq, $counter) = ('', 0);
  while(<IN>){
    chomp;
    my $line = $_;
    $seq.= uc($line) if(eof(IN));
    if (/\>(\S+)/ || eof(IN)){
      if($seq ne ''){
        $counter++;
        if(length($seq) > $lenlimit){
          my $first = substr($seq, 0, $up_iz);
          my $second = substr($seq, -$up_iz);
          my $newseq = $first."NNN".$second;
          print OUT ">$counter\n$newseq\n";
        }
        else{
          print OUT ">$counter\n$seq\n";
        }
      }
      $seq='';
    }else{
      $seq.=uc($line);
    }
  }
  close IN;
  close OUT;
  return $contigfilesub;
}

###FUNCTION TO PARSE THE EVIDENCE FILE, ONLY USED IF TAB FILE IS INSERTED
#Function determines the position of the contigs on the scaffolds, information is used to update the contigs of the TAB file
sub parseEvidenceFile{
   my ($file) = @_;
   my $track_tigs;
   open(IN,$file) || die "Can't open $file -- fatal\n";
   my $scaf = 0;
   my $totalsize= 0;
   while(<IN>){
      chomp;
      if(/^>/){
        $scaf++;
        $totalsize=0;
      }else{
        my ($tig, $size, $links, $gap, $merge) = split(/\|/,$_);
        if($tig ne ""){
          my ($direction, $tig2) = split(/_tig/,$tig);
          $tigOnScafHash->{$tig2}{'begin'} = $totalsize;

          my (undef, $size2) = split(/size/,$size);
          my $end = $totalsize + $size2;
          if($merge ne ""){
            my (undef, $merge2) = split(/merged/,$merge);
            $totalsize = $totalsize + ($size2 - $merge2);
          }elsif($gap ne ""){
            my (undef, $gap2) = split(/gaps/,$gap);
            $gap2 = 1 if($gap2 < 0);
            $totalsize = $totalsize + $size2 + $gap2;
          }else{
            $totalsize = $totalsize + $size2;
          }
          $tigOnScafHash->{$tig2}{'scaf'} = $scaf;
          $tigOnScafHash->{$tig2}{'end'} = $totalsize;
          $tigOnScafHash->{$tig2}{'direction'} = $direction;
        }
      }
   }
}

###FUNCTION TO UPDATE THE ORIGINAL CONTIG FILE INSERTED BY THE USER, SO MULTIPLE TAB FILES OF SEVERAL LIBRARIES CAN BE INSERTED
sub updateContigs{
  my ($file, $update) = @_;

  &printMessage("\n=>".getDate().": Updating contig file\n");

  my ($countContig, $seq, $prevhead) = (0, "", '');
  open(IN,$file) || die "Can't open $file -- fatal\n";
  while(<IN>){
     my $line = $_;
     chomp $line;
     $seq.= $line if(eof(IN));
     if (/\>(\S+)/ || eof(IN)){
       my $head=$1;
       if($prevhead ne ''){
         ++$countContig;
         $tigHash->{$prevhead} = $countContig;
       }
       $prevhead = $head;
       $seq='';
     }else{
        $seq.=$line;
     }
  }
  CounterPrint("                ");
  &FlushFiles();
}

#READ THE CONTIG TO A HASH AND STORE THIS HASH
sub readFileContigHash{
  my ($file) = @_;

  &printMessage("\n=>".getDate().": Reading contig file\n");
    
  my ($contigs, $tig_length);
  my ($countContig, $seq, $prevhead, $step) = (0, "", '', 1000);
  open(IN,$file) || die "Can't open $file -- fatal\n";
  while(<IN>){
     my $line = $_;
     chomp $line;
     $seq.= $line if(eof(IN));
     if (/\>(\S+)/ || eof(IN)){
       my $head=$1;
       if($prevhead ne ''){
         if(++$countContig == $step){
           CounterPrint($countContig);
           $step = $step + 100000;
         }
         $tig_length->{$countContig} = length($seq);
         $contigs->{$countContig}{'name'} = $prevhead;
         $contigs->{$countContig}{'seq'} = $seq;
       }
       $prevhead = $head;
       $seq='';
     }else{
        $seq.=$line;
     }
  }
  CounterPrint("                ");
  &FlushFiles();
  my $contigstore = "$base_name/tmp.$base_name/contigs.stored"; 
  store \%$contigs, "$contigstore";
  undef $contigs;
  return ($contigstore, $tig_length);
}

###FUNCTION THAT FILTERS OUT THE REPEATS BY FINDING CONTIGS THAT HAVE MULTIPLE LINKS WITH OTHER CONTIGS
sub determineRepeats{
  my ($tig_length, $repeathash) = @_;
  my $removeHash;
  #go through each contig
  foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){
    for(my $i = 0; $i < 2; $i++){
      my $dtig = "r" . $tig;
      $dtig = "f" . $tig if($i);
      my $list = $pair->{$dtig};  #get contig pairs from $tig
      my ($seen_it, $matchhash);
      my $ct=0;
      #Go through each contig pair and get the number of links and gapsize
      foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
        my $matchnum = $1 if($match=~/[fr](\w+)/);  
        print TMP "$dtig has $list->{$match}{'links'} links with $match and gap of ".int($list->{$match}{'gaps'}/$list->{$match}{'links'})." bases\n" if($list->{$match}{'links'} >= $min_links);
        if($list->{$match}{'links'} >= $min_links && !defined $seen_it->{$matchnum} && $ct < 2){
          $ct++;
          $matchhash->{$match}{'links'} = $list->{$match}{'links'};
          $matchhash->{$match}{'gaps'} = $list->{$match}{'gaps'};
          $matchhash->{$match}{'ratio'} = $list->{$match}{'gaps'}/$list->{$match}{'links'};
          $seen_it->{$matchnum}++;
      }
      }
      my @arraymatch;
      foreach my $ratiosort (sort {$matchhash->{$a}{'ratio'}<=>$matchhash->{$b}{'ratio'}} keys %$matchhash){
        push @arraymatch, $ratiosort;
      }
      my $repeat = 1;
      my $used;
      my $nummatch = $#arraymatch;
      #only determine if contig is a repeat if it has more than 1 link with other contigs
      if($nummatch > 0){
        my $listmatch = $pair->{$arraymatch[0]};
        #if the top two pairs of $tig have link with each other, establish their link so they are combined in scaffolding stage
        if($listmatch->{$arraymatch[1]}{'links'} >= $min_links){
          $pair = establishLink($dtig, $arraymatch[0], $pair);
          $pair = establishLink($arraymatch[0], $arraymatch[1], $pair);
        }else{ #otherwise, the contig has multiple links and is likely a repeat
          my @linkmatch;
          foreach my $linksort (sort {$matchhash->{$b}{'links'}<=>$matchhash->{$a}{'links'}} keys %$matchhash){
            push @linkmatch, $linksort;
          }
          my ($ratio2, $first, $second) = (0,"","");

          #check for two ratio's  between the two best contig pairs. One is a ratio of the links, other is the number of links per searchspace. 
          #If either one of the two ratio's is above the user-defined (-a) ratio, the original contig is treated as a repeat
          
          #estimate the ratio of the links of the two best contig pairs (ratio 1)
          my $link1 = $matchhash->{$linkmatch[1]}{'links'};
          my $link2 = $matchhash->{$linkmatch[0]}{'links'};
          my $ratio1 = $link1 / $link2;        ## relative ratio of the two most abundant contig pairs
          $ratio1 = sprintf("%.2f", $ratio1);
          $first = $linkmatch[0];
          #estimate the number of links per gap for the two best contig pairs and divide them (ratio 2)
          my $gapPerSpace1 = estimateLinksPerGap($matchhash, $linkmatch[0], $insert_size, $tig_length);
          my $gapPerSpace2 = estimateLinksPerGap($matchhash, $linkmatch[1], $insert_size, $tig_length);
          if($gapPerSpace1 > $gapPerSpace2){
            $second = $linkmatch[0];
            $ratio2 = $gapPerSpace2/$gapPerSpace1;
          }else{
            $second = $linkmatch[1];
            $ratio2 = $gapPerSpace1/$gapPerSpace2;
          }
          my $revdtig = $dtig;
          $revdtig =~ tr/fr/rf/;
          #if one of the two ratio's is above the user-defined (-a) option, contig is a repeat and all links with this contig are removed
          if($ratio2 >= $max_link_ratio || $ratio1 >= $max_link_ratio || $first ne $second){
            foreach my $linksort (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
              my $num = $1 if($linksort=~/[fr](\w+)/);
              $removeHash->{$dtig}{$linksort}++;
              my $revlinksort = $linksort;
              $revlinksort =~ tr/fr/rf/;
              $removeHash->{$revdtig}{$revlinksort}++;
            }
          }
          else{ #otherwise, establish the link between the most likely contig pair
            foreach my $linksort (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
              if($linksort ne $first){
                my $num = $1 if($linksort=~/[fr](\w+)/);
                my $revlinksort = $linksort;
                $revlinksort =~ tr/fr/rf/;
                $removeHash->{$revdtig}{$revlinksort}++;
                $removeHash->{$dtig}{$linksort}++;
              }
            }
          }
        }
      }
    }
  }
  return $removeHash;
}

###FUNCTION TO ESTABLISH A LINK BETWEEN CONTIGS SO THESE CONTIGS ARE PAIRED DURING SCAFFOLDING
sub establishLink{
  my ($tig1, $tig2, $pair) = @_;
  my $list = $pair->{$tig1};
  my $revtig1 = $tig1;
  $revtig1 =~ tr/fr/rf/;
  my $revtig2 = $tig2;
  $revtig2 =~ tr/fr/rf/;
  foreach my $rep (keys %$list){
    if($rep ne $tig2 && $rep ne $revtig2){
      delete $pair->{$tig1}{$rep};
      $rep =~ tr/fr/rf/;
      delete $pair->{$rep}{$revtig1};
    }
  }
  my $list2 = $pair->{$revtig2};
  foreach my $rep2 (keys %$list2){
    if($rep2 ne $tig1 && $rep2 ne $revtig1){
      delete $pair->{$revtig2}{$rep2};
      $rep2 =~ tr/fr/rf/;
      delete $pair->{$rep2}{$tig2};
    }
  }

  return $pair;
}

###DETERMINE THE NUMBER OF LINKS PER GAP, BASED ON INSERT SIZE
sub estimateLinksPerGap{
  my ($linkhash, $tig1, $insert_size, $length_hash) = @_;
  my $t1 = $1 if($tig1=~/[fr](\w+)/);
  my $space = 0;
  my $gap = int($linkhash->{$tig1}{'ratio'});
  $gap = 0 if($linkhash->{$tig1}{'ratio'} < 0);
  if(($length_hash->{$t1}+$gap) >= $insert_size){
    $space = int($insert_size - $gap);
  }else{
    $space =$length_hash->{$t1};
  }
  my $ratio = $linkhash->{$tig1}{'links'}/$space;
  return $ratio;
}

###FUNCTION TO BUILD THE SCAFFOLDS
sub buildScaffolds{
   my ($pair, $tig_length, $verbose, $scaffold, $lib) = @_;
   &printMessage("\n=>".getDate().": Building scaffolds file\n");

   open (SC, ">$scaffold") || die "Can't write to $scaffold -- fatal\n";
   my ($sc_ct, $keyrep, $numrepeat) = (0,0,0);
   my ($repeathash, $seen_start);

   #determine the repeats and remove any link if contig is a repeat
   #if contig has multiple links, but one considered to be the 'best', establish this contig-pair by removing the links with other contigs
   open (TMP, ">$base_name/intermediate_results/$base_name"."_$library.foundlinks.txt") || die "Can't write to $base_name/intermediate_results/$base_name"."_$library.foundlinks.txt -- fatal\n";
   $repeathash = determineRepeats($tig_length, $repeathash);
   close TMP;
   open (REPEAT, ">$base_name/intermediate_results/$base_name"."_$library.repeats.txt") || die "Can't write to $base_name/intermediate_results/$base_name"."_$library.repeats.txt -- fatal\n";
   foreach my $rep (sort keys %$repeathash){
     my $tig = $1 if($rep=~/[fr](\w+)/);;
     my $ls = $repeathash->{$rep};
     my ($num_match,$repline) = (0,"");
     foreach my $rep2 (sort keys %$ls){
       if($pair->{$rep}{$rep2}{'links'} >= $min_links){
         $repline.="\twith $rep2 (links = $pair->{$rep}{$rep2}{'links'})\n";
         $num_match++;
       }
       delete $pair->{$rep}{$rep2};
       delete $pair->{$rep2}{$rep};
     }
     if($num_match > 1){
       $numrepeat++;
       print REPEAT "Contig $rep (size = $tig_length->{$tig}) has $num_match multiple links;\n";
       print REPEAT "$repline\n";
     }
   }
   close REPEAT;
   print SUMFILE "\nREPEATS: \n";
   print SUMFILE "\tNumber of repeated edges = $numrepeat\n$seplines\n";

   #go through each contig and find contig pairs left and right, forming scaffolds
   SEED:
   foreach my $tig (sort {$tig_length->{$b}<=>$tig_length->{$a}} keys %$tig_length){
      my $ftig = "f" . $tig;
      my $rtig = "r" . $tig;

      if(! defined $seen_start->{$tig}){##should prevent re-using a contig as seed if it's already been incorporated into a scaffold
         CounterPrint(++$sc_ct);
         my $chainleft = "";
         my $ori_chainright = $ftig . "Z" . $tig_length->{$tig};
         my $chainright = $ori_chainright;
         my $total = $tig_length->{$tig};
         ($total, $chainright, $seen_start) = &computeLayout("R", $chainright, $ftig, $pair, $tig_length, $total, $seen_start, $tig);
         ($total, $chainleft, $seen_start) = &computeLayout("L", $chainleft, $rtig, $pair, $tig_length, $total, $seen_start, $tig);

         delete $pair->{$ftig};
         delete $pair->{$rtig};
         delete $tig_length->{$tig};
         $seen_start->{$tig}++;
         my $scaffold = $chainleft . $chainright;
         print SC "scaffold" . $sc_ct . ",$total,$scaffold\n";
      }
   }
   CounterPrint("                ");
   close SC;
   &FlushFiles();
}

# links contigs together into a chain - must satisfy user-defined criterions (-k -a)
sub computeLayout{
   my ($ext, $chain, $tig, $pair, $tig_length, $total, $seen_start, $orig_tig_number) = @_;
   my $orig_tig = $tig;
   my $extension = 1;
   EXTENSION:
   while($extension){
      my $tnum = $1 if($tig=~/[fr](\w+)/);
      my $tnumf = "f" . $tnum;
      my $tnumr = "r" . $tnum;
      my $ratio = 0.00;
      if(!defined $seen_start->{$tnum}){ #if already seen in scaffold, do not use it again
        $seen_start->{$tnum}++ if($tnumf ne $orig_tig);
         my $list = $pair->{$tig};
         my $matchhash;
         my ($match1,$link1,$gaps1,$match2,$link2,$gaps2,$cntloop, $countmatches)=("",0,0,"",0,0,0,0);
         my $ct=0;
         LINK:
         foreach my $match (sort {$list->{$b}{'links'}<=>$list->{$a}{'links'}} keys %$list){
            my $matchnum = $1 if($match=~/[fr](\w+)/);
            if($list->{$match}{'links'} >= $min_links && !defined $seen_start->{$matchnum} && $matchnum ne $orig_tig_number && $ct < 2){
              $ct++;
              $matchhash->{$match}{'links'} = $list->{$match}{'links'};
              $matchhash->{$match}{'gaps'} = $list->{$match}{'gaps'};
              $matchhash->{$match}{'ratio'} = $list->{$match}{'gaps'}/$list->{$match}{'links'};
              $countmatches++;
            }else{
              last LINK;
            }
         }
         my $foundlinks = 0;
         if($countmatches > 1){
           my @arraymatch;
           foreach my $ratiosort (sort {$matchhash->{$a}{'ratio'}<=>$matchhash->{$b}{'ratio'}} keys %$matchhash){
             push @arraymatch, $ratiosort;
           }
           my $nummatch = $#arraymatch;
           for(my $i=0; $i <= $nummatch && $foundlinks < 1; $i++){
             my $listmatch = $pair->{$arraymatch[$i]};
              for(my $j=$i+1; $j <= $nummatch && $foundlinks < 1; $j++){
                 my $linkmatch = $listmatch->{$arraymatch[$j]}{'links'};
                 $foundlinks = 1 if(!($linkmatch >= $min_links));
              }
           }
           my $tignum = $1 if($arraymatch[$nummatch]=~/[fr](\w+)/);
           $countmatches=0 if(!$foundlinks && defined $seen_start->{$tignum});
         }if($foundlinks && $countmatches > 1){
             my @linkmatch;
             foreach my $linksort (sort {$matchhash->{$b}{'links'}<=>$matchhash->{$a}{'links'}} keys %$matchhash){
               push @linkmatch, $linksort;
             }
             my $linkhash;
             my $link1 = $matchhash->{$linkmatch[1]}{'links'};
             my $link2 = $matchhash->{$linkmatch[0]}{'links'};
             my $ratio = $link1 / $link2;        ## relative ratio of the two most abundant contig pairs
             $ratio = sprintf("%.2f", $ratio);

             if($ratio <= $max_link_ratio){
               foreach my $mat (keys %$matchhash){
                 delete $matchhash->{$mat} if($mat ne $linkmatch[0]);
               }
               $foundlinks = 0;
               $countmatches = 1;
             }
         }
         if((!$foundlinks) && $countmatches > 0){
           my $nummatch =0;
           my @chainlist;
           my @tiglist;
           foreach my $incl_matches (sort {$matchhash->{$a}{'ratio'}<=>$matchhash->{$b}{'ratio'}} keys %$matchhash){
             if($tig ne $incl_matches){
               $nummatch++;
               my $listmatch = $pair->{$tig};
               my $tempnum = $1 if($incl_matches =~ /[fr](\w+)/);
               my $link2 = $listmatch->{$incl_matches}{'links'};
               my $mean2 = $listmatch->{$incl_matches}{'gaps'}/$link2;

               $seen_start->{$tempnum}++if($nummatch < $countmatches);

               ($chain, $total, $tig) = &getChain($chain, $ext, $link2, $mean2, $incl_matches, $tempnum, $ratio, $tig_length, $total);
               delete $tig_length->{$tempnum};
             }
           }
           $extension = 1;

         }else{
           $extension = 0;
           last EXTENSION;
         }
      }else{
           $extension = 0;
           last EXTENSION;
      }
   }
   return $total, $chain, $seen_start;
}

###function to combine contigs into a scaffold
sub getChain{
  my ($chain, $ext, $link, $mean, $match, $tempnum, $ratio, $tig_length, $total) = @_;
  my $tig = $match;
  if($ext eq "R"){
               $chain .= "k" . $link . "a" . $ratio . "m" . int($mean) . "_" . $match . "z" . $tig_length->{$tempnum};
  }else{
    my $temp_match = "";
    if($match =~ /^r(\d+)/){$temp_match = "f" . $1;}else{$temp_match = "r". $1;}
     $chain = $temp_match . "z" . $tig_length->{$tempnum} . "k" . $link . "a" . $ratio . "m" . int($mean) . "_" . $chain;
 }
 
  $total += $tig_length->{$tempnum};
  return ($chain, $total, $tig);
}


###GET THE DISTANCE BETWEEN TWO PAIRED READS
sub getDistance{

   my ($insert_size, $length_i, $start_i, $start_j) = @_;

   # L  ------  --------- R
   # i    ->        <-    j
   #      ....  ......    insert_span
   #      ============    insert_size




   my $insert_span = ($length_i - $start_i) + $start_j;
   my $gap_or_overlap = $insert_size - $insert_span;
   return $gap_or_overlap;
}

###Pair contigs based on mapping of two reads
sub pairContigs{
  my ($trackA, $trackB, $read_a, $read_b) = @_;
  my ($tig_a, $A_start, $A_end) = split(/\|/, $trackA);
  my ($tig_b, $B_start, $B_end) = split(/\|/, $trackB);
  my ($ori_1,$ori_2) = split(//, $ori);
  if($ori_1 eq "R"){
    my ($tmp_A_start, $tmp_A_end) = ($A_start, $A_end);
    ($A_start, $A_end) = ($tmp_A_end, $tmp_A_start);
  }
  if($ori_2 eq "F"){
    my ($tmp_B_start,$tmp_B_end) = ($B_start,$B_end);
    ($B_start,$B_end) = ($tmp_B_end,$tmp_B_start);
 }
  my $ftig_a = "f" . $tig_a;
  my $ftig_b = "f" . $tig_b;
  my $rtig_a = "r" . $tig_a;
  my $rtig_b = "r" . $tig_b;
  my $A_length = $tig_length->{$tig_a};
  my $B_length = $tig_length->{$tig_b};
  if ($tig_a != $tig_b){####paired reads located on <> contigs
    ####Determine most likely possibility
    if ($A_start < $A_end){
      if ($B_end < $B_start){####-> <- :::  A-> <-B  /  rB -> <- rA
        my $d = &getDistance($insert_size, $A_length, $A_start, $B_start);
        print "A-> <-B  WITH $tig_a -> <- $tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen, Astart,Bstart\n" if($verbose);
        if($d >= $min_allowed){
          $pair->{$ftig_a}{$ftig_b}{'links'}++;
          $pair->{$ftig_a}{$ftig_b}{'gaps'} += $d;
         # $pair->{$ftig_a}{$ftig_b}{'gapdist'}{$d}++;
          $pair->{$rtig_b}{$rtig_a}{'links'}++;
          $pair->{$rtig_b}{$rtig_a}{'gaps'} += $d;
        #  $pair->{$rtig_b}{$rtig_a}{'gapdist'}{$d}++;
          $ct_ok_pairs++;
        }else{                               
          my $err_pair = $ftig_a . "-". $ftig_b;
          $err->{$err_pair}{'links'}++;
          $err->{$err_pair}{'gaps'} += $d;
          $ct_problem_pairs++;
           print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-B  WITH tig#$tig_a -> $d <- tig#$tig_b, A=$A_length nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
        }
      }else{#### -> -> ::: A-> <-rB  / B-> <-rA
        my $rB_start = $B_length - $B_start;
        my $d = &getDistance($insert_size, $A_length, $A_start, $rB_start);
        print "A-> <-rB  WITH $tig_a -> <- r.$tig_b GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Alen,Astart,rBstart\n"if($verbose);
        if($d >= $min_allowed){
          $pair->{$ftig_a}{$rtig_b}{'links'}++;
          $pair->{$ftig_a}{$rtig_b}{'gaps'} += $d;
        #  $pair->{$ftig_a}{$rtig_b}{'gapdist'}{$d}++;
          $pair->{$ftig_b}{$rtig_a}{'links'}++;
          $pair->{$ftig_b}{$rtig_a}{'gaps'} += $d;
        #  $pair->{$ftig_b}{$rtig_a}{'gapdist'}{$d}++;
          $ct_ok_pairs++;
        }else{
          my $err_pair = $ftig_a . "-". $rtig_b;
          $err->{$err_pair}{'links'}++;
          $err->{$err_pair}{'gaps'} += $d;
          $ct_problem_pairs++;
          print PET "Pairs unsatisfied in distance within a contig pair.  A-> <-rB  WITH tig#$tig_a -> $d <- tig#r.$tig_b, A=$A_length  nt (start:$A_start, end:$A_end) B=$B_length nt (start:$B_start, end:$B_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
        }
      }
    }else{
      if ($B_end > $B_start){####<-  -> ::: B-> <-A / rA -> <- rB
        my $d = &getDistance($insert_size, $B_length, $B_start, $A_start);
        print "B-> <-A  WITH $tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,Bstart,Astart\n" if($verbose);
        if($d >= $min_allowed){
          $pair->{$ftig_b}{$ftig_a}{'links'}++;
          $pair->{$ftig_b}{$ftig_a}{'gaps'} += $d;
         # $pair->{$ftig_b}{$ftig_a}{'gapdist'}{$d}++;
          $pair->{$rtig_a}{$rtig_b}{'links'}++;
          $pair->{$rtig_a}{$rtig_b}{'gaps'} += $d;
        #  $pair->{$rtig_a}{$rtig_b}{'gapdist'}{$d}++;
          $ct_ok_pairs++;
        }else{
          my $err_pair = $ftig_b . "-". $ftig_a;
          $err->{$err_pair}{'links'}++;
          $err->{$err_pair}{'gaps'} += $d;
          $ct_problem_pairs++;
          print PET "Pairs unsatisfied in distance within a contig pair.  B-> <-A  WITH tig#$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
        }
      }else{                          ####<- <-  :::  rB-> <-A / rA-> <-B
        my $rB_start = $B_length - $B_start;
        my $d = &getDistance($insert_size, $B_length, $rB_start, $A_start);
        print "rB-> <-A WITH r.$tig_b -> <- $tig_a GAP $d A=$A_length ($A_start-$A_end) B=$B_length ($B_start-$B_end) Blen,rBstart,Astart\n" if($verbose);
        if($d >= $min_allowed){
          $pair->{$rtig_b}{$ftig_a}{'links'}++;
          $pair->{$rtig_b}{$ftig_a}{'gaps'} += $d; 
        #  $pair->{$rtig_b}{$ftig_a}{'gapdist'}{$d}++;
          $pair->{$rtig_a}{$ftig_b}{'links'}++;
          $pair->{$rtig_a}{$ftig_b}{'gaps'} += $d; 
        #  $pair->{$rtig_a}{$ftig_b}{'gapdist'}{$d}++;
          $ct_ok_pairs++;
        }else{
          my $err_pair = $rtig_b . "-". $ftig_a;
          $err->{$err_pair}{'links'}++;
          $err->{$err_pair}{'gaps'} += $d;
          $ct_problem_pairs++;
          print PET "Pairs unsatisfied in distance within a contig pair.  rB-> <-A WITH tig#r.$tig_b -> $d <- tig#$tig_a, B=$B_length nt (start:$B_start, end:$B_end) A=$A_length nt (start:$A_start, end:$A_end) CALCULATED DISTANCE APART: $d < $min_allowed\n";
        }
      }
    }
  }else{###Clone, paired reads located on the same contig -- could be used to investigate misassemblies
      print "Pair ($read_a and $read_b) located on same contig $tig_a ($A_length nt)\n" if ($verbose);
      my $pet_size = 0;

      if ($A_start > $B_start && ($B_start < $B_end) && ($A_start > $A_end)){    # B --> <-- A
        $pet_size = $A_start - $B_start;
        if($pet_size >= $low_iz && $pet_size <= $up_iz){
           $total_for_median++;
           $track_insert->{$pet_size}++;
           $ct_ok_contig++;
        }else{
          print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
          $ct_iz_issues++;
        }
    }elsif($B_start > $A_start && ($B_start > $B_end) && ($A_start < $A_end)){ # A --> <-- B
      $pet_size = $B_start - $A_start;
      if($pet_size >= $low_iz && $pet_size <= $up_iz){
        $total_for_median++;
        $track_insert->{$pet_size}++;
        $ct_ok_contig++;
      }else{
        print PET "Pairs unsatisfied in distance within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end CALCULATED DISTANCE APART: $pet_size\n";
        $ct_iz_issues++;
      }
    }else{
      $ct_illogical++;
      print PET "Pairs unsatisfied in pairing logic within a contig.  Pair ($read_a - $read_b) on contig $tig_a ($A_length nt) Astart:$A_start Aend:$A_end Bstart:$B_start Bend:$B_end\n";
    }
  }
}

###Print read pairing results to the summary file, including estimation of mean and median insert size
sub printResultsPairing{
   print PET "------------- Putative issues with contig pairing - Summary  ----------------\n";
   foreach my $err_pair (sort {$err->{$b}{'links'}<=>$err->{$a}{'links'}} keys %$err){
      my $mean_iz = 0;
      $mean_iz = $err->{$err_pair}{'gaps'} / $err->{$err_pair}{'links'} if ($err->{$err_pair}{'links'});
      print PET "Pair $err_pair has $err->{$err_pair}{'links'} links and mean distance = $mean_iz\n";
   }
   close PET;

   my $satisfied = $ct_ok_pairs + $ct_ok_contig;
   my $unsatisfied = $ct_problem_pairs + $ct_iz_issues + $ct_illogical;
   my $totalPairsUsed_reads = $totalPairsUsed * 2;

   #write distribution file
   open (CSV, ">$distribution") || die "Can't open $distribution for writing -- fatal";
   my ($total_is, $overal_is,$median_ins, $stdev,$record, $sumX,$sumX2) = (0,0,0,0,0,0,0);
   my $median_bin = int($total_for_median/2);

   foreach my $is (sort {$a<=>$b} keys %$track_insert){
     for(my $i=0;$i<$track_insert->{$is};$i++){
       $record++;
       $sumX += $is;
       $sumX2 += ($is * $is);
       $median_ins = $is if($record >= $median_bin && $median_ins == 0);
     }
     $overal_is += ($is * $track_insert->{$is});
      print CSV "$is,$track_insert->{$is}\n";
   }
   my ($mean_ins,$sigma) = (0,0);
   if($sumX > 0 && $record > 0){
     $mean_ins = int($sumX/$record);
     $sigma = sprintf("%.2f",sqrt($sumX2/$record - $mean_ins*$mean_ins));
   }
   close CSV;

   print SUMFILE "READ PAIRS STATS:\n";
   print SUMFILE "\tAssembled pairs: $totalPairsUsed ($totalPairsUsed_reads sequences)\n";
   print SUMFILE "\t\tSatisfied in distance/logic within contigs (i.e. -> <-, distance on target: $insert_size +/$min_allowed): $ct_ok_contig\n";
   print SUMFILE "\t\tUnsatisfied in distance within contigs (i.e. distance out-of-bounds): $ct_iz_issues\n";
   print SUMFILE "\t\tUnsatisfied pairing logic within contigs (i.e. illogical pairing ->->, <-<- or <-->): $ct_illogical\n";
   print SUMFILE "\t\t---\n";
   print SUMFILE "\t\tSatisfied in distance/logic within a given contig pair (pre-scaffold): $ct_ok_pairs\n";
   print SUMFILE "\t\tUnsatisfied in distance within a given contig pair (i.e. calculated distances out-of-bounds): $ct_problem_pairs\n";
   print SUMFILE "\t\t---\n";
   print SUMFILE "\tTotal satisfied: $satisfied\tunsatisfied: $unsatisfied\n\n";
   print SUMFILE "\n\tEstimated insert size statistics (based on $total_for_median pairs): \n";
   print SUMFILE "\t\tMean insert size = $mean_ins\n";
   print SUMFILE "\t\tMedian insert size = $median_ins\n";
  # print SUMFILE "\t\tInsert size deviation = $sigma\n$seplines\n";

   &FlushFiles();
   return $pair;
}


###Parse output of Bowtie to find pairs
sub readBowtie{
   my ($input) = @_;
   my ($counter, $pair_found, $line, $prevline, $prevread) = (0,0, "","","");
   my ($seq1, $seq2, $track1, $track2, $subPairHash);
   my @prevresult;
   open(IN, "$input") || die "Can't open bowtie output -- fatal\n";
   #go through mapping results
   while($line = <IN>){
      $counter++;
      my @result = split(/\t/,$line);
      if($prevread eq $result[0]){
        $pair_found++;
        ($seq1, $track1) = StoreResultsBowtie(@prevresult);
        ($seq2, $track2) = StoreResultsBowtie(@result);
        $subPairHash->{"$track1,$seq1,$track2,$seq2"}++;
      }
      $prevread = $result[0];
      @prevresult = @result;
   }
   close IN;
   
   $totalReads+=$counter;
   $totalPairsFound+=$pair_found;
   return $subPairHash;
}

#retrieve the real position of the reads since the current position is based on the edge of the read. Reverse start and end if read was on <--- strand
sub StoreResultsBowtie{
  my (undef, $strand, $tig, $start, $seq) = @_;
  my ($startval, $endval) = (0,0);
  $tig++;
  if($start > $up_iz && $tig_length->{$tig} > $lenlimit){
    $start = ($tig_length->{$tig} - ($subContigLim - $start));
  }
  if($strand eq "+"){
    $startval = $start;
    $endval = $start + length($seq);
    return $seq, "$tig"."|$startval"."|$endval";
  }
  $startval = $start + length($seq);
  $endval = $start;
  $seq = reverseComplement($seq);
  return $seq, "$tig"."|$startval"."|$endval";
}


###Parse output of BWA/BWAsw
sub readBwa{
  my ($aligninput, $samse) = @_;
  my ($line1, $line2, $pair_found, $ctsingle) = ("","",0, 0,);
  my ($subPairHash);
  system("$aligninput") if($aligninput ne "");
  open(SAM, "$samse") || die "Can't open bwa output -- fatal\n";
  #go through mapping results

  while($line1 = <SAM>){
    next if($line1 =~ /^@/);
    $line2 = <SAM>;
    my @result1 = split("\t", $line1);
    next if($result1[2] eq "*");
    next if(!($result1[1] == 0 || $result1[1] == 16));
    next if($line1 =~ /XA:Z:/);
    $ctsingle++;

    my @result2 = split("\t", $line2);
    next if($result2[2] eq "*");
    next if(!($result2[1] == 0 || $result2[1] == 16));
    next if($line2 =~ /XA:Z:/);
    $ctsingle++;
    $pair_found++;
    my ($seq1, $track1) = StoreResultsBWA(@result1);
    my ($seq2, $track2) = StoreResultsBWA(@result2);
    $subPairHash->{"$track1,$seq1,$track2,$seq2"}++;
  }
  $totalReads+=$ctsingle;
  $totalPairsFound+=$pair_found;
  close SAM;
  return $subPairHash;
}

#retrieve the real position of the reads since the current position is based on the edge of the read. Reverse start and end if read was on <--- strand
sub StoreResultsBWA{
  my @result = @_;
  my ($startval, $endval) = (0,0);
  if(--$result[3] > $up_iz && $tig_length->{$result[2]} > $lenlimit){
    $result[3] = ($tig_length->{$result[2]} - ($subContigLim - $result[3]));
  }
  if ($result[1] & 16) {
    $startval = $result[3] + length($result[9]);
    $endval = $result[3];
    $result[9] = reverseComplement($result[9]);
    return $result[9], "$result[2]"."|$startval"."|$endval";
  }
  $startval = $result[3];
  $endval = $result[3] + length($result[9]);
  return $result[9], "$result[2]"."|$startval"."|$endval";
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

#########END PairingAndScaffolding.pl