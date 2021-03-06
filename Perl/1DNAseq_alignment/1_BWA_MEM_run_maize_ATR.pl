#!/usr/local/bin/perl -w

###The purpose of this script:
###1.To perform DNAseq short reads alignment using publicly available dataset (as along as we know the SRR number)
###2. Because of the limited computing resources we have, we splitted the fastq.gz file into small segments which also allows us to perform parallel computing
###3. Also because of the limited resources, we only save the bam files for the target regions (The regions we want to search for polymorphisms) 
###extract the aligned result at the Maize ATR (Zm00001d01483) gene region to check the polymorphisms

use strict;


my $PE_pl = 'paired_end_trim.pl';

##specify reference genome file
my $genome_ref_4_bwa = '/XXX/Refs/AGPV4/ZmB73_V4_bwa';

##specify output file
my $pre_dir = '/XXX/Align_result/';
mkdir $pre_dir unless -e $pre_dir;
##command line input to control to align which sample
my ($srs_s, $srs_e) = ($ARGV[0], $ARGV[1]);

my $split_lines = 8000000; 
my $split_filter = '--filter='."'".'gzip > $FILE.gz'."'";

##specify the Zm00001d01483 gene target region (It need to be changed if the targetting genome region is different)
my $awk_std1 = '$3==5&&$4>64068899&&$4<64100088';

##The purpose here is to have that single quote for the awk command 
my $awk_std = "'".$awk_std1."'"; 

my $bwa_group = 1;
my ($sra_arrayref, $sra_hashref) = SRA_list_DIR($bwa_group);
my @srs = @$sra_arrayref; 

##specify the header file
my $sam_header_file = '/XXX/Refs/AGPV4/ZmB73_V4_all_sam_header';

for (my $i = $srs_s; $i <= $srs_e; $i ++) { 
	my $sample = $srs[$i];
	my $sample_dir = $pre_dir.$sample.'/'; 
	mkdir $sample_dir unless -e $sample; 
	my $merged_bam = $pre_dir.$sample.'.bam'; ## if the merged file already exists, then we will skip this one
	
	##Get the remote file
	my @remote_files = @{ $$sra_hashref{$sample}};
	my $srs = $remote_files[-1]; 
	my ($pre_srs, $x) = $srs =~ /(SRR\d\d\d)(\d+)/;
	my $remote_ftp_pre = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/'.$pre_srs.'/'; 
	
	my $k; 
	foreach my $run (@remote_files) {
		$k ++;
		my $strain = $sample.'_'.$k;
		my $run_dir = $sample_dir.$run.'/';
		mkdir $run_dir unless -e $run_dir;
		
		##specify a set of file names which will be used to extract from the sra file
		my $ori_run_sra_file = $run_dir.$run.'.sra';
		my $ln_run_sra_file = $sample_dir.$strain.'.sra';
		my $left_file = $sample_dir.$strain.'_1.fastq';
		my $right_file = $sample_dir.$strain.'_2.fastq';
		my $left_gz_file = $left_file.'.gz';
		my $right_gz_file = $right_file.'.gz';
		my $left_chunk_prex = $sample_dir.$strain.'_1_';
		my $right_chunk_prex = $sample_dir.$strain.'_2_';
		my $run_quality_file = $sample_dir.$strain.'_fastq_sample';
		
		my $remote_ftp_file = $remote_ftp_pre.$srs.'/'.$run.'/'.$run.'.sra';
		
		if ($k > 0) {
			##download .sra file. 'wget -q -nH' turn off wget output and disable generation of host-prefixed directory
			system("wget -q -nH $remote_ftp_file -P $run_dir") unless -e $ori_run_sra_file;
			##'ln -s' create a link to TARGET with the name LINK_NAME, -s make symbolic links instead of hard links
			system("ln -s  $ori_run_sra_file $ln_run_sra_file") unless -e $ln_run_sra_file;
			##split .sra file into two separate fastq files (for paired-end sequencing)
			system("fastq-dump -I -split-files --gzip $ln_run_sra_file -O $sample_dir") unless -e $left_file && -e $right_file;
			##Extract about 1000 reads from the 'XXX_1.fastq.gz'file, keep the original file unchanged
			system("gzip -cd $left_gz_file | head -n 4000 > $run_quality_file") unless -e $run_quality_file;
			
			system("rm $ln_run_sra_file") ;
			system("rm $ori_run_sra_file -f");
			
			###gunzip -c t.gz | split -l 1000 -d -a 3 --filter='gzip > $FILE.gz' - x
			##split files by every 8000000 line and generate suffix of length 3
			system("gunzip -c $left_gz_file | split -l $split_lines -d -a 3 $split_filter - $left_chunk_prex");
			system("rm $left_gz_file");
			system("gunzip -c $right_gz_file | split -l $split_lines -d -a 3 $split_filter - $right_chunk_prex");
			system("rm $right_gz_file");
			
		}
		##Get the quality code by calculating the quality score based on the first 1000 reads
		my $quality_code = Quality_Code($run_quality_file);
		
		for (my $i = 0; $i <= 2000; $i ++) {
			my $suffix = sprintf "%03d", $i;
			my $left_chunk_file = $left_chunk_prex.$suffix.'.gz';
			next unless -e $left_chunk_file;
			my $right_chunk_file = $right_chunk_prex.$suffix.'.gz';
			
			##Specify a set of file names which will be used to store intermediate and final files
			my $trim_left_chunk = $left_chunk_file.'.fq_trm';
			my $trim_left_chunk_s = $trim_left_chunk.'_sum';
			my $trim_right_chunk = $right_chunk_file.'.fq_trm';
			my $trim_right_chunk_s = $trim_right_chunk.'_sum';
			
			my $trim_PE_1_chunk = $trim_left_chunk.'.pe';
			my $trim_PE_1_chunk_s = $sample_dir.$strain.'_1_'.$suffix.'.pe_sum';
			
			my $trim_PE_2_chunk = $trim_right_chunk.'.pe';
			my $trim_PE_2_chunk_s = $sample_dir.$strain.'_2_'.$suffix.'.pe_sum';
			
			my $trim_SE_1_chunk = $trim_left_chunk.'.se';
			my $trim_SE_2_chunk = $trim_right_chunk.'.se';
			my $trim_SE_chunk = $sample_dir.$strain.'_0_'.$suffix.'.se';
		
			my $PE_target_bam = $sample_dir.$strain.'_pe_Target.bam_'.$suffix;
			my $SE_target_bam = $sample_dir.$strain.'_se_Target.bam_'.$suffix;

			my $trim_PE_chunk_sam = $sample_dir.$strain.'_pe_all.sam_'.$suffix;
			my $trim_SE_chunk_sam = $sample_dir.$strain.'_se_all.sam_'.$suffix;
			
			my $trim_PE_chunk_uni_sam = $sample_dir.$strain.'_uni_pe.sam_'.$suffix;
			my $trim_SE_chunk_uni_sam = $sample_dir.$strain.'_uni_se.sam_'.$suffix;
		
			my $trim_PE_chunk_uni_bam = $sample_dir.$strain.'_pe_uni.bam_'.$suffix;
			my $trim_SE_chunk_uni_bam = $sample_dir.$strain.'_se_uni.bam_'.$suffix;
			my $chunk_uni_bam = $sample_dir.$strain.'_uni.bam_'.$suffix;
			
			my $chunk_uni_srt_bam = $sample_dir.$strain.'_Target_'.$suffix.'.bam';
			unlink $left_chunk_file if -e  $chunk_uni_srt_bam;
			unlink $right_chunk_file if -e  $chunk_uni_srt_bam;
			next if -e $chunk_uni_srt_bam;
			unlink $trim_PE_chunk_sam if -e $trim_PE_chunk_sam;
			unlink $trim_SE_chunk_sam if -e $trim_SE_chunk_sam;
			
			##trim the fastq file based on the quality score
			if ($quality_code eq 'I') {
				system("Btrim64 -i -q -t $left_chunk_file -Z -l 40 -o $trim_left_chunk -s $trim_left_chunk_s");
				system("rm $left_chunk_file");
				system("Btrim64 -i -q -t $right_chunk_file -Z -l 40 -o $trim_right_chunk -s $trim_right_chunk_s");
				system("rm $right_chunk_file");
				}
				elsif ($quality_code eq 'S') {
					system("Btrim64 -q -t $left_chunk_file -Z -l 40 -o $trim_left_chunk -s $trim_left_chunk_s");
					system("rm $left_chunk_file");
					system("Btrim64 -q -t $right_chunk_file -Z -l 40 -o $trim_right_chunk -s $trim_right_chunk_s");
					system("rm $right_chunk_file");
					}
					
			##/usr/local/bin/paired_end_trim.pl <summary file 1> <summary file 2> <trim output file 1> <trim output file 2>
			##Post-processing paired-end reads after btrim trimming using 'paired_end_trim.pl' program
			system("$PE_pl $trim_left_chunk_s $trim_right_chunk_s $trim_left_chunk $trim_right_chunk");
			system("rm $trim_left_chunk_s $trim_right_chunk_s $trim_left_chunk $trim_right_chunk");
			
			##If the read only in trim left file1 or trim right file2 passed
			system("cat $trim_SE_1_chunk $trim_SE_2_chunk > $trim_SE_chunk");
			system("rm $trim_SE_1_chunk $trim_SE_2_chunk");
			
			##Using the BWA-MEM algorithm to align the paired ends, then extract the target gene region
			system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_PE_1_chunk $trim_PE_2_chunk | awk $awk_std >> $trim_PE_chunk_sam");
			system("cat $sam_header_file $trim_PE_chunk_sam | samtools view - -q 10 -Sb | samtools sort - -o $PE_target_bam")	;
			system("rm $trim_PE_1_chunk $trim_PE_2_chunk $trim_PE_chunk_sam");
			
			system("bwa mem -t 2 -v 0 $genome_ref_4_bwa $trim_SE_chunk | awk $awk_std >> $trim_SE_chunk_sam");
			system("cat $sam_header_file $trim_SE_chunk_sam | samtools view - -q 10 -Sb | samtools sort - -o $SE_target_bam")	;
			system("rm $trim_SE_chunk $trim_SE_chunk_sam");
			
			##Merge the alignments
			system("samtools merge $chunk_uni_bam $PE_target_bam $SE_target_bam");
			system("rm $PE_target_bam $SE_target_bam");
			
			##Sort alignment file
			system("samtools sort $chunk_uni_bam -o $chunk_uni_srt_bam");
			system("rm $chunk_uni_bam");
			
			}
		
		}
		my $target_bams = $sample_dir.'*Target*bam';	
		##Samtools merge all the targeted bam files together, and overwrite the output file if present
		system("samtools merge $merged_bam $target_bams -f");	
		##Index the coordinate-sorted BAM file for fast random access
		system("samtools index $merged_bam");
	}

##Get the SRA file ID for each of the accession 
sub SRA_list_DIR {
	my ($bwa_group) = @_;
	my (%hash, @array);
	open (F, '/XXX/sraRun_by_SRRs_list') || die;
	while (<F>) {
		chomp;
		my @t = split /\t/;
		next if $t[0] eq 'bwa_Group';
		next unless $t[0] == $bwa_group;  
		push @array, $t[1];
		  @{ $hash{$t[1]} } = @t[3..$#t];  
		}
	close F;
	return (\@array, \%hash);	
	}
	
##Check quality score is coded based on the first 1000 reads
sub Quality_Code {
	my ($f) = @_;
	my ($j, $flag, %hash, $code) ;
	open (F, $f) || die;
	while (<F>) {
		chomp;
		my $line = $_;
		$flag = 0 if $line =~ /^\+/;
		$flag ++;
		if ($flag) {
			$j ++;
			my @t = split //, $line;
			foreach my $t(@t) {
				##ord() takes a character and convert it into its ASCII code
				my $n = ord($t);
				$hash{$n} ++;
				}
			}
		last if $j > 1000	;
		}
	close F;	
	if (exists $hash{90}) { $code = 'I'} elsif (exists $hash{48}) {$code = 'S'};
	return $code;		
	}

