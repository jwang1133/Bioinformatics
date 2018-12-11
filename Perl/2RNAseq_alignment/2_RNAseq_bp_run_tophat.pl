#!/usr/local/bin/perl -w

###The purpose of this script:1.To perform RNAseq short reads alignment using publicly available dataset with TopHat(as along as we know the SRR number)
###1. This script is designed for paired-end sequencing reads
###2. Because of the limited resources, we only save the bam files for the target regions (The regions we want to search for polymorphisms) 
###3. We only aligned result at the Gmax ligase1 gene region due to the limited resources

use strict;
##specify the master directory
my $dir = '/XXX/RNASeq/';
my $srr_list_file = $dir.'SRR_list';

##infer the SRR code by calling function Parse_srr_list()
my $srr_arrayref = Parse_srr_list($srr_list_file);
my $remote_sra_pre = 'ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByStudy/sra/SRP/SRP048/SRP048326/';

##specify reference file
my $ref_4_bowtie2 = '/XXX/Refs/Gmax_RefGen_V2/Gmax_RefGen_V2';

##loop through each of the SRR to perform the alignment
foreach my $srr (@$srr_arrayref) {

	my $srr_dir = $dir.$srr.'/';
	mkdir $srr_dir unless -e $srr_dir;
	my $srr_tophat_dir = $srr_dir.'tophat/';
	mkdir $srr_tophat_dir unless -e $srr_tophat_dir;
	my $remote_file = $remote_sra_pre.$srr.'/'.$srr.'.sra';
	my $local_file = $srr_tophat_dir.$srr.'.sra';
	
	##download files using wget
	system("wget  -q $remote_file  -O $local_file");
	
	##split .sra file into two separate fastq files (for paired-end sequencing reads)
	system("fastq-dump -I -split-files --gzip $local_file -O $srr_tophat_dir");
	system ("rm $local_file");
	
	my $srr_1_file = $srr_tophat_dir.$srr.'_1.fastq.gz';
	my $srr_2_file = $srr_tophat_dir.$srr.'_2.fastq.gz';
	my $accepted_bam = $srr_tophat_dir.'accepted_hits.bam';
	my $lig1_bam = $srr_dir.'Lig1.bam';
	
	##Perform alignment with tophat2
	system("tophat2 -p 10 -o $srr_tophat_dir $ref_4_bowtie2 $srr_1_file $srr_2_file");
	system("samtools index $accepted_bam");
	
	##Extract the target bam file region with samtools
	system("samtools view $accepted_bam 11:26629471-26638425 -q 50 -b > $lig1_bam");
	system("rm -r $srr_tophat_dir");
	
	}


#######Function to pass the SRR code from the file
sub Parse_srr_list {
	my ($f) = @_;
	open (F, $f) || die;
	my @array;
	while (<F>) {
		chomp;
		my $line = $_;
		push @array, $line;
		}
	close F;
	return \@array;	
	}