#!/usr/bin/perl

# open virus genome FASTA
open(VF,"/datastore/nextgenout4/share/labs/bioinformatics/Selitsky/simulated_human_and_viral_reads/viral_genome_unmasked.fa");

# store in a hash
$count = 0;
while($header = <VF>){
   chomp $header;
   $seq = <VF>;
   chomp $seq;
   $count = $count + 1;
   $header =~ s/>/\@/g;
   $info = join("\t",$header,$seq);
   $virHashCount{$count} = $info;
}

@simLen = (1000);
@mut = (0,1,2,3,4,5,6,7,8,9,10);
for($j=1;$j<101;$j++){
   for($k = 0;$k<scalar(@simLen);$k++){
      for($i=0;$i<$simLen[$k];$i++){
	 $ranVir = int(rand(1893));
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

### add random mut
	    if($read1 !~ /N/ && $read2 !~ /N/){
	       $base{0} = "A";
	       $base{1} = "T";
	       $base{2} = "C";
	       $base{3} = "G";
	       for($n=0;$n<scalar(@mut);$n++){
		  @seq = split("",$read1);
		  %mutLocHash = ();
		  for($m=0;$m<$mut[$n];$m++){
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
		  open(OUT1,">>/datastore/nextgenout4/share/labs/bioinformatics/Selitsky/simulated_human_and_viral_reads/mutVirSim/all_vir_numSim_".$simLen[$k]."_iter_".$j."_numMut_".$mut[$n]."_R1.fastq");
		  open(OUT2,">>/datastore/nextgenout4/share/labs/bioinformatics/Selitsky/simulated_human_and_viral_reads/mutVirSim/all_vir_numSim_".$simLen[$k]."_iter_".$j."_numMut_".$mut[$n]."_R2.fastq");
		  $header =~ s/\//_/g;
		  print OUT1 "$header:$ranSeq\n$newSeq\n+\nQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n";
		  print OUT2 "$header:$ranSeq\n$read2\n+\nQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ\n";
	       }
	    }
	    else{
	       $i=$i-1;
	    }
	 }
	 else{
	    $i=$i-1;
	 }
      }
   }
}
