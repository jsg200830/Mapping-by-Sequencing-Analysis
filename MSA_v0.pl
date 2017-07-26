#!/usr/bin/perl

# Copyright (c) 2015
# University of Nebraska - Lincoln
# All Rights Reserved
#
# Authors: Shangang Jia, David Holding and Chi Zhang
#
# Redistribution and use in source and binary forms, with or without modification, are permitted provided that the 
# following conditions are met:
#   o Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
#   o Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following 
#     disclaimer in the documentation and/or other materials provided with the distribution.
#   o Neither the name of the University of Nebraska - Lincoln nor the names of its contributors may be used to endorse or promote 
#     products derived from this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED 
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A 
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SHANGANG JIA, DAVID HOLDING and CHI ZHANG 
# (OR UNIVERSITY OF NEBRASKA - LINCOLN) BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR 
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, 
# STRICT LIABILITY, OR  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS 
# SOFTWARE, EVEN IF ADVISED OF THE  POSSIBILITY OF SUCH DAMAGE.
#
# DESCRIPTION: MSA enables a complete analysis of mapping-by-sequencing, by calling the thrid party softwares, including Samtools,
# Bowtie2 and VarScan. It is an integrated pipeline from input of raw fastq data, alignment, SNPs/indel calling, and plotting mapping linkage peaks.
# It supports inputs of fastq, bam, or SNP files in a comparison of mutant type and its contrast normal type
# for F2 generation. This script generates the results in one command line.
#
# CITATION: it can be cited as Shangang Jia, David Holding and Chi Zhang, unpublished.
#


use warnings;
use Getopt::Long;
use Time::Local;
use Cwd;
use constant DEFAULT_TEMP_DIRECTORY => ".";				# Default directory
my $pwd = cwd();
# Generate a run ID to be used for temporary files
sub generateID {
    my ($sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst) = localtime(time);
    return sprintf("%d%02d%02d.%02d%02d%02d", 1900+$year, $mon+1, $mday, $hour, $min, $sec);
} # End of sub generateID

# Given total seconds signifying the module run-time, format it in hh:mm:ss format
sub formatTime {
    my $totaltime = $_[0];
    my $str;

    $str = sprintf("%02d:", $totaltime / 3600);    # total hours
    $totaltime = $totaltime % 3600;
    $str .= sprintf("%02d:", $totaltime / 60);      # total minutes
    $totaltime = $totaltime % 60; 
    $str .= sprintf("%02d", $totaltime);            # total sconds

    return $str;
} # End of sub formatTime

# CNS windows moving on the chromosomes, with a windows of 100kb, step of 100kb.
sub densitymoving {
    my ($bamdir, $filenamex, $mutantnum, $wildnum, $myreference) = @_;
    my ($win, $step);
    if (!$win){$win=100000;}
    if (!$step){$step=100000;}
    
    # retrieve chromosome length from the reference.
    my (%chrlength,@f1,@fields);
    system "bowtie2-inspect -s $myreference > $bamdir/ref.summary.txt";
    unless ( open (PRO, "$bamdir/ref.summary.txt") ) {
        print "Cannot open file \"$bamdir/ref.summary.txt\"\n";
        exit;
    }
    
    <PRO>;
    while(<PRO>) {
        chomp($_);
        next if ($_ !~ /Sequence/);
        @fields = split('\t',$_);
        @f1 = split(' ',$fields[1]);
        $chrlength{$f1[0]}=$fields[2];
    }
    close PRO;
    # retrieve chromosome length from the reference.
    
    open (OUTPUTa, ">$bamdir/mu_wt.posall");
    print OUTPUTa "chr,start,positivecns,allcns,ratio\n";
    open (OUTPUTb, ">$bamdir/mu_wt.pos");
    print OUTPUTb "chr,start,ref,var,samples_as_provided\n";
    print "Started, processing \"$filenamex, win-step, $win,$step\"\n";
    #reading all cns frequencies.
    unless ( open (PROx, $filenamex) ) {
        print "Cannot open file \"$filenamex\"\n";
        exit;
    }
    <PROx>;
    my (%allcnsx,$posx,@fields);
    while(<PROx>) {
        chomp($_);
        @fields = split(',',$_);
        $posx = $fields[0]."_".$fields[1];
        $allcnsx{$posx}=$_;
    }
    close PROx;
    my ($chlekey, $positivenum,$negativenum,$iallnum,$key,$ratio,$pos,$position,$istart,$allsamples,$e,$j,$k,$l,@mfre,@wfre,$msflag,$mbflag,$wsflag,$wbflag,$fma,$fwa);
    ###for each chromosome.

    for $chlekey (sort keys %chrlength) {
    $position = 0;
    $istart=0;
        #go from the start to the end of one chromosome.
        until ($position>$chrlength{$chlekey}){
        $position = $istart+$win;
        $positivenum=0;
        $negativenum=0;
        $iallnum=0;
            for ($e=$istart; $e<= $position; $e++) {
            $key = $chlekey."_".$e;
            
            ####if cns exists.
                if ($allcnsx{$key}){
                @fields = split(',',$allcnsx{$key});
                shift @fields; shift @fields; shift @fields; shift @fields;
                    $allsamples=$mutantnum+$wildnum;
                    for ($j=0; $j<$#fields; $j+=2) {
                        $k = int($j/2);
                        $l =$j+1;
                        if ($k<$mutantnum){$fma = $fields[$j]/($fields[$j]+$fields[$l]);push @mfre, $fma;}
                        if ($k>=$mutantnum){$fwa = $fields[$j]/($fields[$j]+$fields[$l]);push @wfre, $fwa;}
                    }
                    $msflag=1;$mbflag=1;
                    for ($j=0; $j<=$#mfre; $j++) {
                        if ($mfre[$j]<0.05){$mbflag=0;}
                        elsif ($mfre[$j]>0.95){$msflag=0;}
                        else {$mbflag=0;$msflag=0;}
                    }
                    $wsflag=1;$wbflag=1;
                    for ($j=0; $j<=$#wfre; $j++) {
                        if ($wfre[$j]<0.5){$wbflag=0;}
                        elsif ($wfre[$j]>=0.5){$wsflag=0;}
                        else {$wbflag=0;$wsflag=0;}
                    }
                    undef @mfre;undef @wfre;
                    
                    if ($msflag && $wbflag){
                        $positivenum ++;$pos .= $allcnsx{$key}."\n";
                    }
                    elsif ($mbflag && $wsflag){
                        $positivenum ++;$pos .= $allcnsx{$key}."\n";
                    }
                    else {$negativenum ++;}
                $iallnum ++;	
                }
            ####if cns exists.
            
            }
        if (!$positivenum){$positivenum=0;}
        if (!$negativenum){$negativenum=0;}
        if (!$iallnum){$iallnum=0;}
        $ratio = $positivenum/($iallnum+0.05);
        print OUTPUTa "$chlekey,$istart,$positivenum,$iallnum,$ratio\n";
            print "";
        if ($pos){print OUTPUTb "$pos";}
        undef $pos;
        $istart += $step; 
        }
        #go from the start to the end of one chromosome.
    }
   ###for each chromosome.
    close OUTPUTa;
    close OUTPUTb;
    
}
#end of densitymoving.

# chromosome plot by ggbio.
sub chromplotggbio {
	my ($bamdir, $plotfilename, %chrlength) = @_;
    my ($firstch,$mychname,$mychle);
    for my $chlekey (sort keys %chrlength) {
        $firstch ="chr".$chlekey;
        $mychname =$chlekey;
        $mychle = $chrlength{$chlekey};
    }
	open (OUTPUT, ">$bamdir/ggbio.R");
	print OUTPUT "";
    
	print OUTPUT "
	library(ggbio)
	library(GenomicRanges)

	gr <- GRanges(seqnames =\"$firstch\", IRanges(start = 1, end = 2))
	seqinfo(gr) <- Seqinfo(paste0(\"chr\", c($mychname)), c($mychle), NA, \"mock1\")
	xq <- seqinfo(gr)
	p1<-autoplot(xq[paste0(\"chr\", c($mychname))])

	setwd(\"$bamdir/\")
	d=read.csv(\"$plotfilename\", header=TRUE)
	start=d[,2]
	end=d[,2]+100000
	chr=d[,1]
	den=d[,3]
	snp <- GRanges(seqnames =paste0(\"chr\", data.matrix(as.numeric(as.character(chr)))), IRanges(start = data.matrix(as.numeric(as.character(start))), end = data.matrix(as.numeric(as.character(end)))), den=den)
	p2<-p1+layout_karyogram(snp, aes(x = start, y = den), geom = \"line\", ylim = c(0, 10), color = \"red\", fill =\"red\")
	pdf(\"$plotfilename.ggbio.pdf\")
	p2
	dev.off()
	
	";
	close OUTPUT;

	system "Rscript $bamdir/ggbio.R";	
}

# Samtool for mileup file, which is processed by VarScan for cns files with fitlering.
sub mileup2varscan2cns {
	my ($bamdir, $reference, $allmutantbam, $allwildbam) = @_;
	#mileup by samtools.
    my ($j, $mutantnum,$wildnum,$allsamnum,@mbam,@wbam);
    print "$reference\n$allmutantbam\n$allwildbam\n!!!!!!!!!\n\n\n";
	system "samtools mpileup -f $reference $allmutantbam $allwildbam > $bamdir/mutant_wild.mileup";
    @mbam = split(' ',$allmutantbam);
    @wbam = split(' ',$allwildbam);
    $allsamnum =scalar(@mbam)+scalar(@wbam);
	#SNP/InDel calling by VarScan.	
	print STDERR sprintf("\n\n");
	print STDERR sprintf("Progressing: VarScan beggins\n");
	print STDERR sprintf("\n");

	system "java -jar VarScan.v2.3.7.jar mpileup2cns $bamdir/mutant_wild.mileup > $bamdir/mutant_wild.mileup.cns";
	my ($snppo, $refvar, @fre, $fre, @gspli, $sumsingle, $refacc, $varacc);
	my $wt_mu_cns = "$bamdir/mutant_wild.mileup.cns";
    #save the filtered cns file.
    open (OUTPUT, ">$bamdir/mutant_wild.mileup.cns.fil");
    print OUTPUT "\n";
    
	unless ( open (WTMUCNS, $wt_mu_cns) ) {
		print "Cannot open file \"$wt_mu_cns\"\n";
		exit;
	}
	while (<WTMUCNS>) {
   
        chomp;
		@gspli = split("\t",$_);
		$snppo = $gspli[0].",".$gspli[1];#cns position with chr and position.
		if ($gspli[2]=~ m/,/){next;}#discard the alternative type for reference.
		if ($gspli[3]=~ m/,/){next;}#discard the alternative type for variation.
		$refvar = $gspli[2].",".$gspli[3];#ref and var.
		@fre = split(":",$gspli[10]);#all frequency data in 11th column.
		#set unfound position's fre to 0.
        $refacc=2;
        $varacc=3;
        if ($fre[$refacc] eq "-"){$fre[$refacc]=0;}
        if ($fre[$varacc] eq "-"){$fre[$varacc]=0;}
        $fre = $fre[$refacc].",".$fre[$varacc];#frequency data for wt and mu.
        if ($fre[$refacc]=~ /^[+-]?\d+$/){}else {
            print STDERR sprintf("\n\n");
            print STDERR sprintf("WARNING: Not an integer for ref allele frequecy. $fre[$refacc]\n$gspli[10]\n");
            print STDERR sprintf("\n");
            next;
        }
        if ($fre[$varacc]=~ /^[+-]?\d+$/){}else {
            print STDERR sprintf("\n\n");
            print STDERR sprintf("WARNING: Not an integer for var allele frequecy. \n");
            print STDERR sprintf("\n");
            next;
        }
        $sumsingle =$fre[$refacc]+$fre[$varacc];
        if ($sumsingle==0){goto mynextfre;}
          for ($j=2; $j<= $allsamnum; $j++) {
              $refacc+=5;$varacc+=5;
              if ($fre[$refacc] eq "-"){$fre[$refacc]=0;}
              if ($fre[$varacc] eq "-"){$fre[$varacc]=0;}
              if ($fre[$refacc]=~ /^[+-]?\d+$/){}else {
                  print STDERR sprintf("\n\n");
                  print STDERR sprintf("WARNING: Not an integer for ref allele frequecy. \n");
                  print STDERR sprintf("\n");
                  next;
              }
              if ($fre[$varacc]=~ /^[+-]?\d+$/){}else {
                  print STDERR sprintf("\n\n");
                  print STDERR sprintf("WARNING: Not an integer for var allele frequecy. \n");
                  print STDERR sprintf("\n");
                  next;
              }
              $sumsingle =$fre[$refacc]+$fre[$varacc];
              if ($sumsingle==0){goto mynextfre;}
              $fre .= ",".$fre[$refacc].",".$fre[$varacc];
          }

        print OUTPUT "$snppo,$refvar,$fre\n";
    mynextfre:
	}
	close OUTPUT;
	close WTMUCNS;
        print "CNS calling finished\n\n\n\n";
}
#end of mileup2varscan2cns.

# define command line arguments
my (@m_fqp, @w_fqp, @m_bam, @w_bam, @m_cns, @w_cns, $reference, $reference2, $output, $cpus);

my $result = &GetOptions("m_fqp:s{0,}" => \@m_fqp,
						 "w_fqp:s{0,}" => \@w_fqp,
						 "m_bam:s{0,}" => \@m_bam,
						 "w_bam:s{0,}" => \@w_bam,
						 "mw_cns:s{0,}" => \@mw_cns,
						 "reference=s{1}" => \$reference,
						 "reference2=s{1}" => \$reference2,
						 "output|o=s{1}" => \$output,
						 "threads|t=i{0,}" => \$cpus,
						 "log!" => \$log,);
## Print help info.
unless ($result && ( (scalar(@m_bam)+scalar(@w_bam))>0 || (scalar(@m_fqp)+scalar(@w_fqp))>0 || scalar(@mw_cns)>0 ) && (defined($reference) || defined($reference2)) && defined($output)) {
	print STDERR sprintf("\n");
	print STDERR sprintf("MSA commands\n");
	print STDERR sprintf("=============================================================\n");
	print STDERR sprintf("USAGE:\n");
	print STDERR sprintf("  perl %s [--m_bam <mutant bam files>] [--w_bam <wild type bam files>] [--m_fqp <mutant fastq paired files>] [--w_fqp <wild type fastq paired files>]\n", $0);
	print STDERR sprintf("[--mw_cns <mutant cns files>] [--threads|-t <CPU threads>] --reference <genome reference fasta file> --output|-o <output.prefix> [OPTIONS]\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("WHERE:\n");
	print STDERR sprintf("  --m_bam <mutant bam files>         : Path to alignment output files generated by bowtie2. The specified files should be only for mutants, which are seperated by space.\n");
	print STDERR sprintf("  --w_bam <wild type bam files>      : Path to alignment output files generated by bowtie2. The specified files should be only for wild type, which are seperated by space.\n");
	print STDERR sprintf("  --m_fqp <fastq paired files for mutant>     : Path to the fastq files with pair-end reads for mutants.\n");
	print STDERR sprintf("  --w_fqp <fastq paired files for wild type>  : Path to the fastq files with pair-end reads for wild type.\n");
    print STDERR sprintf("  The above fq options are only required if mapping by bowtie2 needed, unless bam files are provided.\n");
	print STDERR sprintf("  --mw_cns <mutant cns files>        : Path to the cns file harboring mutant and wild.\n");
	print STDERR sprintf("  --reference <genome reference fasta file>    : Path to the fasta reference file.\n");
    print STDERR sprintf("  --reference2 <bowtie 2 indexed database>    : Path to the bowtie 2 indexed database.\n");
	print STDERR sprintf("  --threads|-t <CPU cores or threads>   : Multiple threads\n");
	print STDERR sprintf("  --output|-o <output.prefix>           : Path to output files' prefix\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("OPTIONS:\n");
	print STDERR sprintf("  --log|--nolog                     : Enable/Disable the cleaning of temporary files [DEFAULT: --log]\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("Example: perl MSA.pl --w_fqp=wt_1.fastq wt_2.fastq --m_fqp=mu_1.fastq mu_2.fastq --reference2=reference.bw2 --threads=20 --output=test \n");
	print STDERR sprintf("BOWTIE2 commands:\n");
	print STDERR sprintf("  %% bowtie2 -x reference -1 fastq1 -2 fastq2 -S output.sam\n");
	print STDERR sprintf("\n");
	print STDERR sprintf("\n");
	exit();
} # End of unless statement

# assigning default values
$cpus = 1 if (!defined($cpus) || $cpus !~ m/^\d+$/);
$log = "true" if (!defined($log));

my $m_mapping_input = scalar(@m_bam);
my $w_mapping_input = scalar(@w_bam);
my $mw_mapping = $m_mapping_input + $w_mapping_input;
my $m_fastq_input = scalar(@m_fqp);
my $w_fastq_input = scalar(@w_fqp);

my $m_cns_input = scalar(@m_cns);
my $w_cns_input = scalar(@w_cns);

#check if there are paired input files.
if ($m_fastq_input == 0 && $w_fastq_input == 0) {
    print STDERR sprintf("\n\n");
    print STDERR sprintf("Warning: No fastq input file specified, you must specify bam or cns input files without fastq input files. \n");
    print STDERR sprintf("\n");
}
elsif (($m_fastq_input == 0 && $w_fastq_input > 0) || ($m_fastq_input > 0 && $w_fastq_input == 0)){
    print STDERR sprintf("\n\n");
    print STDERR sprintf("Warning: No paired fastq input file specified. \n");
    print STDERR sprintf("\n");
}
# End of if statmenet

if ($m_mapping_input == 0 && $w_mapping_input == 0) {
	print STDERR sprintf("\n\n");
	print STDERR sprintf("Warning: No bam input file specified, you must specify fastq or cns input files without bam input files. \n");
	print STDERR sprintf("\n");
}
elsif (($m_mapping_input == 0 && $w_mapping_input > 0) || ($m_mapping_input > 0 && $w_mapping_input == 0)){
    print STDERR sprintf("\n\n");
    print STDERR sprintf("Warning: No paired bam input file specified. \n");
    print STDERR sprintf("\n");
}
# End of if statmenet

if ($m_cns_input == 0 && $w_cns_input == 0) {
	print STDERR sprintf("\n\n");
	print STDERR sprintf("Warning: No cns input file specified, you must specify fastq or bam input files without cns input files. \n");
	print STDERR sprintf("\n");
}
elsif (($m_cns_input == 0 && $w_cns_input > 0) || ($m_cns_input > 0 && $w_cns_input == 0)){
    print STDERR sprintf("\n\n");
    print STDERR sprintf("Warning: No paired cns input file specified. \n");
    print STDERR sprintf("\n");
}
# End of if statmenet



# generate running ID and running folders.
my $run_id = &generateID();
system "mkdir $run_id";
system "mkdir $run_id/bamsam";
my $resultdir = "$pwd/$run_id";
my $bamdir = "$pwd/$run_id/bamsam";
my ($e_ffile_wt, $e_ffile_mu, $b_ffile_wt, $b_ffile_mu);#bam files for merging use in samtools.
my $starttime;
$starttime = timelocal(localtime(time));	#local current time.

#Indexing bowtie2 database.
my ($myreference);

if (defined($reference)) {
    print STDERR sprintf("\n\n");
    print STDERR sprintf("Processing: Indexing bowtie2 database\n");
    print STDERR sprintf("\n");
    system "bowtie2-build $reference $bamdir/reference.bw2 1>$bamdir/log.txt";
    $reference = "$pwd/$reference";
    $myreference = "$bamdir/reference.bw2";
}
elsif (defined($reference2)) {
    print STDERR sprintf("\n\n");
    print STDERR sprintf("Working with provided index file\n");
    print STDERR sprintf("\n");
    system "bowtie2-inspect $reference2 > $bamdir/reference.fa";
    $reference = "$bamdir/reference.fa";
    $myreference = $reference2;

}
else {
    print STDERR sprintf("\n\n");
    print STDERR sprintf("No reference files provided!!!\n");
    print STDERR sprintf("\n");
    exit;
}
#Indexing bowtie2 database.

#check if bam files provided and mapping needed.
my ($j, $k, $n, @msortedbamfiles, @wsortedbamfiles,$sbamname,@mbamfiles,@wbamfiles);
if ($m_mapping_input > 0 && $w_mapping_input > 0) {
    print STDERR sprintf("\n\n");
    print STDERR sprintf("WARNING: bam files provided, and skipping bowtie2 step. \n");
    print STDERR sprintf("\n");
    for ($j=0; $j<= $#m_bam; $j++) {
        push @mbamfiles, $m_bam[$j];
    }
    
    for ($j=0; $j<= $#w_bam; $j++) {
        push @wbamfiles, $w_bam[$j];
    }
    goto mybam;
} # End of if statmenet
if ($m_cns_input > 0 && $w_cns_input > 0) {
    print STDERR sprintf("\n\n");
    print STDERR sprintf("WARNING: cns files provided, and skipping bowtie2 and cns calling steps. \n");
    print STDERR sprintf("\n");

    system "cp $mw_cns[0] $bamdir/mutant_wild.mileup.cns.fil";

    goto mycns;
} # End of if statmenet




#running mapping.
#pair end alignment for mutant.

for ($j=0; $j< $#m_fqp; $j++) {
    $k=$j*2+1;
    $n=$j*2;
    system "bowtie2 -p $cpus -x $myreference -1 $m_fqp[$n] -2 $m_fqp[$k] -S $bamdir/mutant.PE$j.sam 1>$bamdir/log.txt";
    system "samtools view -bS -@ $cpus $bamdir/mutant.PE$j.sam -o $bamdir/mutant.PE$j.bam  1>$bamdir/log.txt";
    $sbamname="$bamdir/mutant.PE$j.bam";
    push @mbamfiles, $sbamname;
}
#pair end alignment for wild.

for ($j=0; $j< $#w_fqp; $j++) {
    $k=$j*2+1;
    $n=$j*2;
    system "bowtie2 -p $cpus -x $myreference -1 $w_fqp[$n] -2 $w_fqp[$k] -S $bamdir/wild.PE$j.sam 1>$bamdir/log.txt";
    system "samtools view -bS -@ $cpus $bamdir/wild.PE$j.sam -o $bamdir/wild.PE$j.bam  1>$bamdir/log.txt";
    $sbamname="$bamdir/wild.PE$j.bam";
    push @wbamfiles, $sbamname;
}
mybam:

for ($j=0; $j<= $#mbamfiles; $j++) {
    system "samtools sort -@ $cpus $mbamfiles[$j] $bamdir/mutant.PE$j.sorted";
    system "samtools index -@ $cpus $bamdir/mutant.PE$j.sorted.bam";
    $sbamname="mutant.PE$j.sorted.bam";
    push @msortedbamfiles, $sbamname;
}
for ($j=0; $j<= $#wbamfiles; $j++) {
    system "samtools sort -@ $cpus $wbamfiles[$j] $bamdir/wild.PE$j.sorted";
    system "samtools index -@ $cpus $bamdir/wild.PE$j.sorted.bam";
    $sbamname="wild.PE$j.sorted.bam";
    push @wsortedbamfiles, $sbamname;
}
#end of mapping for wt and mu.

# retrieve chromosome length from the reference.
my (%chrlength);
print "$myreference\n";
system "bowtie2-inspect -s $myreference > $bamdir/ref.summary.txt";
unless ( open (PRO, "$bamdir/ref.summary.txt") ) {
    print "Cannot open file \"$bamdir/ref.summary.txt\"\n";
    exit;
}

<PRO>;
while(<PRO>) {
    chomp($_);
    next if ($_ !~ /Sequence/);
    @fields = split('\t',$_);
    @f1 = split(' ',$fields[1]);
    $chrlength{$f1[0]}=$fields[2];
}
close PRO;
# retrieve chromosome length from the reference.

# call for SNP/InDel (CNS).
my ($allmutantbam, $allwildbam);
for ($j=0; $j<= $#msortedbamfiles; $j++) {
    $allmutantbam .= "$bamdir/$msortedbamfiles[$j]"." ";
}
for ($j=0; $j<= $#wsortedbamfiles; $j++) {
    $allwildbam .= "$bamdir/$wsortedbamfiles[$j]"." ";
}
my $mutantnum=scalar(@msortedbamfiles);
my $wildnum =scalar(@wsortedbamfiles);
    &mileup2varscan2cns($bamdir, $reference, $allmutantbam, $allwildbam);

# cns density moving.
mycns:
my $filteredcns = "$bamdir/mutant_wild.mileup.cns.fil";
	&densitymoving($bamdir, $filteredcns, $mutantnum, $wildnum, $myreference);

# plotting in the genome coordinates.
my $plotfilename = "mu_wt.posall";
	&chromplotggbio($bamdir, $plotfilename, %chrlength);


my $total_seconds;
	my $endtime = timelocal(localtime(time));
	my $diff = $endtime - $starttime;
	$total_seconds += $diff;		# Accumulate total run-time
	print STDERR sprintf("\n\n");
	print STDERR sprintf("Parsing CNS finished, used time: %s\n", &formatTime($diff));	
	print STDERR sprintf("\n");
	
$endtime = timelocal(localtime(time));
$diff = $endtime - $starttime;
print STDERR sprintf("\n\n");
print STDERR sprintf("All finished, total time: %s\n", &formatTime($diff));	
print STDERR sprintf("\n");	
	
