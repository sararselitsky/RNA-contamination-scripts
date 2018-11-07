#!/usr/bin/perl

# read in hg38 RSEM transcripts FATA
open(VF,"/home/selitsky/ref/hg38/hg38_rsem.transcripts.fa");


$base{0} = "T";
$base{1} = "C";
$base{2} = "G";
$base{3} = "A";

@lowComplexReads = ();
@read=();

#### generate low complexity reads ####

for($i=0;$i<4;$i++){
   for($k=0;$k<50;$k++){
      $read[$k] = $base{$i};
   }
   $newSeq = join("",@read[0..49]);
   push @lowComplexReads, $newSeq;
}

for($i=0;$i<4;$i++){
   for($j=0;$j<4;$j++){
      if($base{$j} ne $base{$i}){
	 for($k=0;$k<50;$k++){
	    $read[$k] = $base{$i};
	    $read[$k+1] = $base{$j};
	    $k=$k+1;
	 }
	 $newSeq = join("",@read[0..49]);
	 push @lowComplexReads, $newSeq;
      }
   }
}

for($i=0;$i<4;$i++){
   for($j=0;$j<4;$j++){
      for($n=0;$n<4;$n++){
	 if($base{$j} ne $base{$i} || $base{$i} ne $base{$n} || $base{$j} ne $base{$n}){
	    for($k=0;$k<50;$k++){
	       $read[$k] = $base{$i};
	       $read[$k+1] = $base{$j};
	       $read[$k+2] = $base{$n};
	       $k=$k+2;
	    }
	    $newSeq = join("",@read[0..49]);
	    push @lowComplexReads, $newSeq;
	 }
      }
   }
}

### store FASTA into hash ####
#
$count = 0;
while($header = <VF>){
   chomp $header;
   $seq = <VF>;
   chomp $seq;
   $seq = uc($seq);
   $count = $count + 1;
   $header =~ s/>/\@/g;
   $info = join("\t",$header,$seq);
   $virHashCount{$count} = $info;
}

### simulate reads #####
@mut = (0,1,2,3,4,5,6,7,8,9,10);
for($j=1;$j<101;$j++){
   open(OUT1,">/datastore/nextgenout4/share/labs/bioinformatics/Selitsky/simulated_human_and_viral_reads/humanSim/humanSim_".$j."_R1.fastq");
   open(OUT2,">/datastore/nextgenout4/share/labs/bioinformatics/Selitsky/simulated_human_and_viral_reads/humanSim/humanSim_".$j."_R2.fastq");
   for($i=0;$i<1000000;$i++){
      $ranVir = int(rand(181937));
      my($header,$seq) = split("\t",$virHashCount{$ranVir});
      if(length($seq) > 600){
	 $ranSeq = int(rand(length($seq)));
	 if($ranSeq > 300){
	    $ranSeq = $ranSeq - 301;
	 }
	 $read1 = substr($seq,$ranSeq,50);
	 $ranSeqPair2 = $ranSeq + 200;
	 $read2 = substr($seq,$ranSeqPair2,50);
	 $read2 = reverse $read2;
	 $read2 =~ tr/ATCG/TAGC/;
	 if($read1 !~ /N/ && $read2 !~ /N/){
	    $header =~ s/\//_/g;
	    print OUT1 "$header:$ranVir:$ranSeq\n$read1\n+\nQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n";
	    print OUT2 "$header:$ranVir:$ranSeq\n$read2\n+\nQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n";
	 }
	 else{
	    $i=$i-1;
	 }
      }
      else{
	 $i=$i-1;
      }
   }

   for($n=0;$n<1000;$n++){
      $randomSeq = int(rand(scalar(@lowComplexReads)));      
      @seq = split("",$lowComplexReads[$randomSeq]);
      %mutLocHash = ();
      $numMut = int(rand(scalar(@mut)));
      for($m=0;$m<$mut[$numMut];$m++){
	 $mutLoc = int(rand(49));
	 $mutBase = int(rand(3));
	 if($seq[$mutLoc] ne $base{$mutBase} && !exists $mutLocHash{$mutLoc}){
	    $seq[$mutLoc] = $base{$mutBase};
	    $mutLocHash{$mutLoc} = 1;
	 }
	 else{
	    $m = $m-1;
	 }
      }
      $newSeq = join("",@seq[0..49]);
      $read2 = reverse $lowComplexReads[$randomSeq];
      $read2 =~ tr/ATCG/TAGC/;
      print OUT1 "\@low_complex:$n:$mut[$numMut]\n$newSeq\n+\nQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n";
      print OUT2 "\@low_complex:$n:$mut[$numMut]\n$read2\n+\nQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n";
   }
}
