# input data format
#7778   	1      	110919_UNC13-SN749_0114_BB007YABXX_8_ACAGTG    	gi|9626069|ref|NC_001357.1|Humanpapillomavirus-18,completegenome       	C

# input files into an array
@files = glob("/datastore/nextgenout4/share/labs/bioinformatics/seqware/kirc_virus/virus-output/*/*/viralRTable.txt");

# loop through files
for($i=0;$i<scalar(@files);$i++){
   open(IN,$files[$i]);
   while($line = <IN>){
      chomp $line;
      my($pos,$cov,$sample,$virus,$nt) = split(/\t/,$line);
      if($nt ne "." && $cov > 1){
# HeLa SNV table
	 $posHash{104} = "C";
	 $posHash{287} = "G";
	 $posHash{485} = "C";
	 $posHash{549} = "A";
	 $posHash{751} = "T";
	 $posHash{806} = "A";
	 $posHash{1012} = "T";
	 $posHash{1194} = "A";
	 $posHash{1353} = "A";
	 $posHash{1807} = "C";
	 $posHash{1843} = "G";
	 $posHash{2269} = "T";
	 $posHash{5875} = "A";
	 $posHash{6401} = "G";
	 $posHash{6460} = "G";
	 $posHash{6625} = "G";
	 $posHash{6842} = "G";
	 $posHash{7258} = "A";
	 $posHash{7486} = "T";
	 $posHash{7529} = "A";
	 $posHash{7567} = "C";
	 $posHash{7592} = "C";
	 $posHash{7670} = "T";

# determine major alleles
	 @nts = split("",$nt);
	 %ntHash = ();
	 for($j=0;$j<$cov;$j++){
	    $letter = $nts[$j];
	    $ntHash{$letter} += 1;
	 }
	 $undef = $major;
	 foreach $key (keys %ntHash){
	    if($ntHash{$key} > $cov/2 ){
	       $major = $key;
	    }
	 }

	 if(exists $posHash{$pos}){
	    $countHash{$sample} += 1;
	    if($major eq $posHash{$pos}){
	       $HelA_Hash{$sample} += 1; 
	    }
	    if($major ne $posHash{$pos}){
	       $not_HelA_Hash{$sample} += 1; 
	    }
	 }
      }
   }
}

# output results
print "sample\ttotal_HPV18_SNPs\tHeLa_SNPS\tref_SNPs\n";
foreach $key (keys %countHash){
   print "$key\t$countHash{$key}\t$HelA_Hash{$key}\t$not_HelA_Hash{$key}\n";
}


