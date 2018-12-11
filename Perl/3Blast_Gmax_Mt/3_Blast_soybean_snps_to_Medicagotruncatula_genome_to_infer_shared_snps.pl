#!/usr/local/bin/perl -w

###The purpose of this script is to find the shared SNPs between soybean amd Medicago truncatula (M.t.) by blast soybean snps to the Medicago truncatula genome
###The script is designed to perform parallel computing with the limited resources

use strict;
use Bio::DB::Fasta;

##specify a set of directories for the input file
my $ref_dir ='/XXX/Refs/Gmax_v2/';
my $ori_dir = '/XXX/genotype_in_segment/';

##speficy the M.t. database file
my $Mt_db = '/XXX/Mtruncatula_285_Mt4.0.fa';

##specify directory for intermediate and final output files
my $odir1 = '/XXX/snp_flanking_seq/';
mkdir $odir1 unless -e $odir1;
my $odir2 = '/XXX/Results/';
mkdir $odir2 unless -e $odir2;

##parameter that will be used for the 
my $parameter = "\"m D\"";

##Command line inputs which will be used for parallel computing
my ($ch_s, $ch_e, $seg_s, $seg_e) = ($ARGV[0], $ARGV[1], $ARGV[2], $ARGV[3]);

##Refer the length of each chromosome by calling function Parse_chro_length()
my $at_chro_size_hashref = Parse_chro_length();

##loop through each of the 20 chromosome files
for (my $ch = $ch_s; $ch <= $ch_e; $ch ++) {
	my $ch_fas_file = $ref_dir.'Gm'.$ch.'.fa';
	
	##using chromosome fasta file as database with Bio::DB::Fasta
	my $ch_fas_db = Bio::DB::Fasta->new($ch_fas_file);
	my $chro_size = $$at_chro_size_hashref{'Chr'.$ch};	
	
	##loop through the genotype file in segments for each chromosome
	for (my $seg = $seg_s; $seg <= $seg_e; $seg = $seg+2 ){
		my $ori_seg_file = $ori_dir.'Gm302_chr'.$ch.'_'.$seg.'Mb';
		next unless -e $ori_seg_file;
		my $blast_result = $odir1.'Gm302_chr'.$ch.'_'.$seg.'Mb.snp_blast_result';
		my $seg_ancestral_allele_file = $odir2.'Gm302_Chr'.$ch.'_'.$seg.'Mb.ancestral_allele';
		open (F, $ori_seg_file) || die;
		open (OUTPUT, '>'.$seg_ancestral_allele_file) || die;
		while (<F>) {
			chomp;
			my @t = split /\t/;		
			next if $t[0] =~ /Chromosome/;	
			my $site = $t[1];	
			my $output_file = $odir1.'Gm302_chr'.$ch.'_'.$seg.'Mb.flanking_seq';
			##specify the start and end position of soybean DNA sequence we want to extract
			my ($site_p, $site_a) = ($site - 29, $site + 29);	
			next unless (($site_p >= 0) && ($site_a <= $chro_size));
			##extract the seqeunce
			my $ref_seq_fas = $ch_fas_db -> seq($ch, $site_p => $site_a);
			
			##Format the extracted soybean DNA sequence file into a blast input file format 
			my $snp_flanking_region_file = Parse_flanking_seq_file($ch, $site, $ref_seq_fas, $output_file);
			##perform blast with blastn
			system("blastall -p blastn -d  $Mt_db -i $snp_flanking_region_file -m 4 -e 0.1 -o $blast_result -F $parameter -b 10");
			
			##Extract the Mt allele from the blast result file.
			open (IN, $blast_result) || die;	 
			my ($line_num,$pline_n) = (0, 1);
			my ($flag, $flag_b) = (1, 0);
			my ($non_base_len, $all_count_len, $snp_position, $Mt_allele);
			while(<IN>){
				chomp;
				my $line = $_; 
				$line_num ++;
				next unless ($line_num >=17);
				if ($line =~ /No hits found/){
					$flag = 0;	
				};
				
				next unless $flag != 0;	
				if ($line =~ /1_0/){
					$pline_n = $line_num;
					$flag ++;
					my $count = 0;
					my @t1 = split //, $line;
					my @t2 = split /\s+/, $line;
					
					##To make sure the matched region including the SNP position 
					if (($t2[1] <= 30) && ($t2[3] >= 30)){
						$flag_b = 1;
						my $flag_a = 1;
						
						##To calculate the position of the SNP we want to extract
						foreach my $e (@t1){
							my $eb = $e; 				
							if($eb =~ /[acgtn]/){
								$flag_a = 0;
							}
							next unless $flag_a != 0;
							if ($eb  !~ /[acgtn]/){
								$count++;
							}
						}
						$non_base_len = $count;
						my @qseq = split//, $t2[2];
						my $all_count = 0;
						my $b_count = 0;
						foreach my $b (@qseq){
							next unless $b_count < 30 - $t2[1] + 1;
							if ($b =~ /[acgtn]/){
								$b_count ++;
							}
							$all_count ++;
						}
						$all_count_len = $all_count;
						$snp_position = $non_base_len + $all_count_len;	 
					}
				}	
				##To make sure only extract Mt allele from the top hit blast result
				next unless (($flag == 2)&& ($flag_b == 1)); 
				if ($line_num - $pline_n == 1){
					##calculate the coverage in bases
					my @t3 = split /\s+/, $line;
					my $coverage = abs($t3[3] - $t3[1]) + 1;
					##extract the the Mt allele 
					$Mt_allele = substr($line,$snp_position-1,1);
					next unless $Mt_allele =~ /[acgt]/;
					$Mt_allele = uc $Mt_allele;
					print OUTPUT 	$ch."\t".$site."\t".$Mt_allele."\t".$coverage."\n";
				}	
			}
			close IN;
		}
		close F;
 }
}

##Function to pass the length of each chromsome
sub Parse_chro_length {
	my $file = '/XXX/Gm302_chromosome_length_base'; 
	my %hash;
	open (FILE, $file) || die;
	while (<FILE>) {
		chomp;
		my @t = split /\t/;
		$hash{$t[0]} = $t[1];
		}
	close FILE;
	return \%hash;
	}	

##Function to format the extracted soybean DNA sequence file into a blast input file format
sub Parse_flanking_seq_file{
	my ($chromosome, $bp, $flanking_sequence, $out_put_file) = @_;
	open (OUT, '>'.$out_put_file) || die;
	print OUT ">chr".$chromosome.'snp'.$bp."\n";
 	print OUT $flanking_sequence."\n";
 	close OUT;
 	return $out_put_file;
	}
	
	
	
