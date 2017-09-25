#!/usr/bin/perl

# exome-seq analysis, relating depth to exons in gtf file. 

use warnings;



###########reading gtf file.
$filenamex = "Zea_mays.AGPv3.25.gtf";
unless ( open (PRO, $filenamex) ) {
print "Cannot open file \"$filenamex\"\n\n";
exit;
}
@cnvr = <PRO>;
close PRO;

for ($we=0; $we<= $#cnvr; $we++) {
		@fs = split ("\t", $cnvr[$we]);
		$chr=$fs[0];
		$exon=$fs[2];
		$start=$fs[3];
		$end=$fs[4];
		$incgene=$fs[8];
		@forid = split('gene_id \"',$incgene);
		@forid2 = split('\"',$forid[1]);
		$id = $forid2[0];
		@foren = split('exon_number \"',$incgene);
		@foren2 = split('\"',$foren[1]);
		$en = $foren2[0];
		if ($exon eq "exon"){
			$postart = $chr."_".$start;
			$exonspos{$postart}=$end;
			$exonstart{$postart}=$id."_E".$en;
		}
	
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

###########reading depth file.
$filenamex = $ARGV[0];
open (OUTPUT, ">$filenamex.exondepth");
print OUTPUT "";

unless ( open (PRO, $filenamex) ) {
print "Cannot open file \"$filenamex\"\n";
exit;
}

<PRO>;
$currchr = 1;
$currsite = 1;
while(<PRO>) {

    chomp($_);
    @fields = split('\t',$_);
	$temv = $fields[2];
		next if ($fields[0] =~ /\D/);
##### same chromosome
        if ($fields[0] eq $currchr){
			for ($i=$currsite; $i<= $fields[1]; $i++) {
				
				$pos = $fields[0]."_".$i;
				if ($exonstart{$pos}){undef $exonaccudepth;undef $currexonpos;undef $currdep;$currexonpos = $pos; $currdep = $exonstart{$pos}.",".$fields[0].",".$i;}
				#####
				if ($currexonpos && $exonspos{$currexonpos} && $i<$exonspos{$currexonpos}){
					if ($temv){$exonaccudepth += $temv;undef $temv;}
				}
				elsif ($currexonpos && $exonspos{$currexonpos} && $i==$exonspos{$currexonpos}) {
					if ($temv){$exonaccudepth += $temv;undef $temv;}
					@forlg = split(',',$currdep);
					$lg = $exonspos{$currexonpos}-$forlg[2];
					if (!$exonaccudepth){$exonaccudepth=0;}
					if ($lg<=0){$avgv=0;print "0\t$currdep,$exonspos{$currexonpos}\t$exonaccudepth\n";}else{
					$avgv = $exonaccudepth/$lg;}
					print OUTPUT "$currdep,$exonspos{$currexonpos},$lg,$exonaccudepth,$avgv\n";
					undef $exonaccudepth; undef $currexonpos;undef $currdep;
				}
				#####
			}
$currsite = $fields[1]+1;
		}
##### same chromosome

##### changed chromosome
		elsif ($fields[0] ne $currchr){$currsite = 1;$pos = $fields[0]."_".$fields[1];
	        if ($exonstart{$pos}){undef $exonaccudepth;undef $currexonpos;undef $currdep;$currexonpos = $pos; $currdep = $exonstart{$pos}.",".$fields[0].",".$fields[1]; }
			goto mynext;
		}
##### changed chromosome	

mynext:
$currchr = $fields[0];
}
close PRO;


close OUTPUT;


