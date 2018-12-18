## This folder describes the process for performing RNAseq short reads alignment on publicly available RNAseq short reads dataset

### I. Things need to know about the script '2\_RNAseq\_bp\_run\_tophat.pl'

**A. This script uses publicly available RNAseq data, directly download from website (as long as we know the SRR code )**
**B. This script is designed to perform alignment for paired-end sequencing reads**
**C. We only save the bam files at the targeted regions at the final step because of the limited resources**

### II. The process for the RNAseq alignment with TopHat

1. Prerequisite tools bowtie2, TopHat, samtools, fastq-dump need to be installed 

	a, 	bowtie2

		Bowtie 2 version 2.2.8 by Ben Langmead (langmea@cs.jhu.edu, www.cs.jhu.edu/~langmea)

	b, 	samtools
		
		Program: samtools (Tools for alignments in the SAM format)
		Version: 1.6 (using htslib 1.6)


	d,  fastq-dump

		fastq-dump : 2.8.2

2. Prepare references with bowtie2-build 

	a, Download the reference genome fasta file (As we only interested in the chromosome regions, so download only chromosome.*.fa.gz)

		$ cd /XXX/Refs/Gmax_RefGen_V2/
		$ for i in {1..20}; do  wget  ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.0.dna.chromosome."$i".fa.gz ; done 
		$ for i in {1..20}; do cat Glycine_max.Glycine_max_v2.0.dna.chromosome."$i".fa.gz  >> Gmax_RefGen_V2.fa.gz ; done
		$ gunzip Gmax_RefGen_V2.fa.gz

	b, Build the references with bowtie2-build

		$ bowtie2-build /XXX/Refs/Gmax_RefGen_V2/Gmax_RefGen_V2.fa /XXX/Refs/Gmax_RefGen_V2/Gmax_RefGen_V2

3. Prepare SRR list, example SRR_list looks like this:

		SRR1592289
		SRR1592290
		SRR1592291
		SRR1592292
		SRR1592293

4. Run the perl script named '2\_RNAseq\_bp\_run\_tophat.pl'

		$perl /XXX/scripts/2_RNAseq_bp_run_tophat.pl &

		