# Mapping-by-Sequencing-Analysis
Bulk segregation analysis mapping is a widely used approach in functional genomics, especially now that low sequencing costs make mapping-by-sequencing realistic. However, the huge data processing required is not trivial. To help in this regard, we utilized our population of maize kernel mutants to develop a new mapping-by-sequencing analysis software tool, (MSA), to identify causative genes based on the DNA and RNA level variations in mutant pools. It is designed to identify the association regions, which show linkage disequilibrium to causative mutations, and place linkage peaks on the chromosomes. The method involves the comparison of mutant and normal pools of F2 individuals, and the parents for F2s can be also used alongside the two pools. It allows simultaneous analysis of multiple mutants and multiple normal samples, which can reduce the noise. 

# Dependencies required
1.	Perl models, Getopt::Long, Time::Local, and Cwd. 
2.	Samtools. 
3.	Bowtie2.
4.	Java and VarScan. 
5.	R base and Bioconductor package “ggbio”.

# Input
The data necessary to run MSA consists of the reference genome and three types of input data, fastq sequencing data, bam files, and SNP/indel files, for one of the following comparisons: 
1. mutant v.s. normal pool, where each pool has at least 20 individuals. The normal pool may be pure wild type or a mixture of wild type homozygotes and heterozygotes. 
2. add parents for more reliable positive SNPs which reduce the noise to make linkage peaks. For example, (mutant pool & B73) vs (wild type pool & Mo17) in F2s of B73 mutants x Mo17. 
3. multiple mutant pools vs multiple wild type pools, for example, introducing more wild type SNPs to reduce the noise and make linkage peak more accurate. 

The example input files are listed as followed:

Fastq paired input files (recommended, MSA will run bowtie2 alignment firstly, and produce final output):

mu_1.fastq & mu_2.fastq: pair end fastq files for mutant pool.

wt_1.fastq & wt_2.fastq: pair end fastq files for wild type pool.

Bam paired input files (recommended, MSA will not run alignment step but go directly into SNP/indel refining for postive ones, and produce final output):

mutant.PE0.sorted.bam: aligment file for mutant pool by bowtie2.

wild.PE0.sorted.bam: aligment file for wild type pool by bowtie2.

Variant input file (time effective and only for linkage peak plotting, but need to organize SNP/indel file):

mutant_wild.mileup.cns.fil: SNP/indel file with a format of chromosome#, position, reference nucleotide, alternative nucleotide, reference read number, alternative site' read number. It is retrieved from the SNP/indels file by Varscan, and the user can make it by themselves by using any available software. 

# Running examples
Start with paired fastq files (wild type VS mutant type)

1.	“perl MSA_v0.pl --w_fqp=wt_1.fastq wt_2.fastq --m_fqp=mu_1.fastq mu_2.fastq --reference2=reference.bw2 --threads=20 --output=test”.

Start with alignemnt bam files which should be paired for wild type and mutant type.

2.	“perl MSA_v0.pl --m_bam=mutant.PE0.sorted.bam --w_bam=wild.PE0.sorted.bam --reference2=reference.bw2 --threads=20 --output=test”.

Start with a SNP/indel file which should be organized in a format as followed: chromosome#, position, reference nucleotide, alternative nucleotide, reference read number, alternative site' read number (see example input file "mutant_wild.mileup.cns.fil").

3.	“perl MSA_v0.pl --mw_cns=mutant_wild.mileup.cns.fil --threads=20 --output=test”

# Output
MSA displays the linkage peaks, which indicate linkage disequilibrium and causative mutations responsible for the phenotype of interest. Results consist of SNPs/indels, positive SNPs/indels, and linkage peak figures. It allows simultaneous analysis of multiple mutants, multiple wild type pools, and does the comparison between them, by searching for the identical regions in the mutants against wild type. For example, in F2s from B73 mutant x Mo17 maize mapping populations, the mutants would be screened for the B73 SNPs/indels, which are considered as positive SNPs. The Mo17 SNPs/indels in the mutant are discarded, which reduces the noise in the whole genome. To extend the search based on more samples, any mutant pools may be used as an input in MSA, compared to the wild type without causative mutations and mutant phenotype. 

The example output files are listed as followed:

mu_wt.pos: all SNPs/indels produce by Varscan. Format: chromsome#, position, reference nucleotide, alternative nucleotide, reference read number for wild type, alternative site' read number for wild type, reference read number for mutant, alternative site' read number for mutant

mu_wt.posall: all postive SNPs/indels which are homozygous in the mutant pool and heterzygous in the wild type pool, in a window of 100kb. Format: chromsome#, position, positive SNPs/indels, all SNPs/indels, frequency of postive SNPs/indels. 

mu_wt.posall.ggbio.pdf: linkage peaks plotted by R. 

# Summary for all files in example folder
Fastq files: Wild type, wt_1.fastq and wt_2.fastq; Mutant type, mu_1.fastq and mu_2.fastq.

Bam files: wild.PE0.sorted.bam and mutant.PE0.sorted.bam.

SNP file: mutant_wild.mileup.cns.fil.

Reference: reference.fa and bowtie2 index files: reference.bw2.*. 

Result files: mu_wt.posall.ggbio.pdf (linkage peak plot), mu_wt.pos (postive SNPs/indels), mu_wt.posall (positve SNP number in a window of 100kb)

# Deletion identification
In addition, we provided two perl scripts for deleted exon identification: exome_depth_maize.pl and comparing_depth_maize.pl. They can work with bam files from wild type and mutant type samples, and produced a table to compare the abundance of exons based on the B73 assembly (for example v3.25, downloaded from ftp://ftp.ensemblgenomes.org/pub/plants/release-25/gtf/zea_mays/Zea_mays.AGPv3.25.gtf.gz).

Runnig steps:
1. Alignment files by bowtie2, similar to the above MSA. (wild.PE0.sorted.bam and mutant.PE0.sorted.bam)
2. download maize gtf file, ftp://ftp.ensemblgenomes.org/pub/plants/release-25/gtf/zea_mays/Zea_mays.AGPv3.25.gtf.gz.
3. gunzip Zea_mays.AGPv3.25.gtf.gz
4. samtools depth wild.PE0.sorted.bam > wild.PE0.sorted.bam.depth
5. samtools depth mutant.PE0.sorted.bam > mutant.PE0.sorted.bam.depth
6. perl exome_depth_maize.pl wild.PE0.sorted.bam.depth
7. perl exome_depth_maize.pl mutant.PE0.sorted.bam.depth
8. perl comparing_depth_maize.pl wild.PE0.sorted.bam.depth mutant.PE0.sorted.bam.depth

The output file, wild.PE0.sorted.bam.depth.exondepth_mutant.PE0.sorted.bam.depth.exondepth_exondepth.txt, is with read depth for each exon with wild type vs mutant type as followed.
ExonID  Chromosome# Start End Length  Depth_wildtype  Normalizeddepth_wildtype  ExonID  Chromosome# Start End Length  Depth_mutanttype  Normalizeddepth_mutanttype  annotation

The gene annotation is from the file "biomart_maize_genefeature.txt". 




