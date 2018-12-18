## This folder describes the process for performing DNAseq short reads alignment on publicly available DNAseq short reads dataset

### I. Things need to know about the script '1\_BWA\_MEM\_run\_maize\_ATR.pl'

**A. This script uses publicly available DNAseq data, directly download from website**

**B. We split the fastq.gz file into small segments because of the limited computing resources, this also allow us to perform parallel computing**

**C. We only save the bam files at the targeted regions at the final step because of the limited resources**

### II. The process for the DNAseq alignment

1. Prerequisite tools BWA, samtools, Btrim64, fastq-dump need to be installed 

	a, 	BWA

		Program: bwa (alignment via Burrows-Wheeler transformation)
		Version: 0.7.5a-r405
		Contact: Heng Li <lh3@sanger.ac.uk>

	b, 	samtools
		
		Program: samtools (Tools for alignments in the SAM format)
		Version: 1.6 (using htslib 1.6)

	c, 	Btrim64

		Author: Yong Kong

		Reference: Kong, Y (2011) Btrim: A fast, lightweight adapter and quality trimming program for next-generation sequencing
		 technologies, Genomics, 98, 152-153.

		Contact: yong.kong@yale.edu

	d,  fastq-dump

		fastq-dump : 2.8.2


2. Build reference genome index using bwa index (Linux command line)

		$wget -o /XXX/Refs/AGPV4/Zea_mays.AGPv4.dna.toplevel.fa.gz   ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz &
		$gunzip /XXX/Refs/AGPV4/Zea_mays.AGPv4.dna.toplevel.fa.gz
		$cd /XXX/Refs/AGPV4/
		$ln -s Zea_mays.AGPv4.dna.toplevel.fa ZmB73_V4_bwa
		$bwa index -a bwtsw ZmB73_V4_bwa

3. Build index for fasta file

		$mv /XXX/Refs/AGPV4/Zea_mays.AGPv4.dna.toplevel.fa /XXX/Refs/AGPV4/ZmB73_V4_all.fa
		$samtools faidx /XXX/Refs/AGPV4/ZmB73_V4_all.fa

4. Create sequence dictionary file

		$samtools dict ZmB73_V4_all.fa -o ZmB73_V4_all_sam.dict

		Extract the sam header file

		$cut -f1-3 ZmB73_V4_all_sam.dict > ZmB73_V4_all_sam_header

5. Prepare the SRR list, example 'sraRun\_by\_SRRs\_list' looks like following：

		bwa_Group	Sample_Name_s	Experiment_s	SRRs
		1	B73-1	SRX131321	SRR447984
		1	B97-1	SRX131294	SRR447957
		1	BKN015	SRX131088	SRR447751
		1	BKN031	SRX131193	SRR447856
		1	CML69	SRX131156	SRR447819
		1	TIL01	SRX131219	SRR447882
		1	TIL05	SRX131092	SRR447755
		2	BKN009	SRX131279	SRR447942
		2	BKN025	SRX131158	SRR447821
		2	BKN026	SRX131199	SRR447862
		2	TIL07	SRX131297	SRR447960
		2	TIL11	SRX131246	SRR447909

6. Run the perl script named '1\_BWA\_MEM\_run\_maize\_ATR.pl' 

		$perl /XXX/scripts/1_BWA_MEM_run_maize_ATR.pl 0 4 &