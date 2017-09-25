#!/usr/bin/perl

# comparing the exon depth in exome-seq analysis. 

use warnings;



###########reading gene feature file.
$filenamex = "biomart_maize_genefeature.txt";
unless ( open (PRO, $filenamex) ) {
print "Cannot open file \"$filenamex\"\n\n";
exit;
}
@cnvr = <PRO>;
close PRO;
shift @cnvr;
for ($we=0; $we<= $#cnvr; $we++) {
		chomp $cnvr[$we];
		@fs = split (",", $cnvr[$we]);
			$genef{$fs[0]}=$fs[1];
	
}



$filenamex = $ARGV[0];
@forsa = split('_',$filenamex);
unless ( open (PRO, $filenamex) ) {
print "Cannot open file \"$filenamex\"\n\n";
exit;
}
@cnvr = <PRO>;
close PRO;

for ($we=0; $we<= $#cnvr; $we++) {
		chomp $cnvr[$we];
		@fs = split (",", $cnvr[$we]);
		$key = $fs[1]."_".$fs[2];
			$exon1{$key}=$cnvr[$we];
	
}

$filenamex = $ARGV[1];
@forsa2 = split('_',$filenamex);
unless ( open (PRO, $filenamex) ) {
print "Cannot open file \"$filenamex\"\n\n";
exit;
}
@cnvr = <PRO>;
close PRO;

for ($we=0; $we<= $#cnvr; $we++) {
	chomp $cnvr[$we];
		@fs = split (",", $cnvr[$we]);
		$key = $fs[1]."_".$fs[2];
			$exon2{$key}=$cnvr[$we];
	
}

for $key (keys %exon1)
    {
		@fs = split (",", $exon1{$key});
		if (!$exon2{$key}){$exon2{$key}=$fs[0].",".$fs[1].",".$fs[2].",".$fs[3].","."0".","."0".","."0";}

	}
	
for $key (keys %exon2)
	    {
			@fs = split (",", $exon2{$key});
			if (!$exon1{$key}){$exon1{$key}=$fs[0].",".$fs[1].",".$fs[2].",".$fs[3].","."0".","."0".","."0";}

		}
	
		$chr1 = 301433382;
		$chr2 = 237893627;
		$chr3 = 232227970;
		$chr4 = 242029974;
		$chr5 = 217928451;
		$chr6 = 169381756;
		$chr7 = 176810253;
		$chr8 = 175347686;
		$chr9 = 157021084;
		$chr10 = 149627545;
		$chr1 = 301433382;
		$chr2 = 237893627;
		$chr3 = 232227970;
		$chr4 = 242029974;
		$chr5 = 217928451;
		$chr6 = 169381756;
		$chr7 = 176810253;
		$chr8 = 175347686;
		$chr9 = 157021084;
		$chr10 = 149627545;	
		open (OUTPUT, ">$forsa[0]\_$forsa2[0]_exondepth.txt");
		print OUTPUT "";		
		for ($w=1; $w<= 10; $w++) {
			$forch = "chr".$w;
			print "$$forch\n";
				for ($e=1; $e<= $$forch; $e++) {
					$key = $w."_".$e;
					####if exon exists.
					if ($exon1{$key}){
						@geneid = split(',',$exon1{$key});
						@geneid2 = split('_',$geneid[0]);
						$geneid = $geneid2[0];
						if (!$genef{$geneid}){$genef{$geneid} = "NA";}
					print OUTPUT "$exon1{$key},$exon2{$key},$genef{$geneid}\n";
					}		
					####if exon exists.
		
			   }

		    }


		close OUTPUT;
