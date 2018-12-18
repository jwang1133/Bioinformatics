## This folder describes the process for finding shared SNPs between soybean amd Medicago truncatula (M.t.) by blast soybean snps to the Medicago truncatula genome

### I. Things need to know about the script '3\_Blast\_soybean\_snps\_to\_Medicagotruncatula\_genome\_to\_infer\_shared\_snps.pl'

**A. The purpose of this script is to find the shared SNPs between soybean amd Medicago truncatula (M.t.) by blast soybean snps to the Medicago truncatula genome**

**B. The script is designed to perform parallel computing with the limited resources**

### II. The process for the blast

1. Prerequisite tools BLASTN (Version:2.2.28) need to be installed

2. Prepare the Mtuncatula database file

	a, Download the Mtuncatula database file from 'http://genome.jgi.doe.gov/pages/dynamicOrganismDownload.jsf?organism=Phytozome'

	b, Decompress the file
		
		$ gunzip Mtruncatula_285_Mt4.0.fa.gz


	c, Format the database 

		$ formatdb -p F -i Mtruncatula_285_Mt4.0.fa -o T

3. Download the Gmax reference genome fasta file

		$ cd /XXX/Refs/Gmax_RefGen_V2/
		$ for i in {1..20}; do  wget  ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/glycine_max/dna/Glycine_max.Glycine_max_v2.0.dna.chromosome."$i".fa.gz ; done 
		$ gunzip Glycine_max.Glycine_max_v2.0.dna.chromosome.*.fa.gz
		$ for i in {1..20}; do mv Glycine_max.Glycine_max_v2.0.dna.chromosome."$i".fa  Gm"$i".fa; done &

4. Genotype file has been processed as in .hmp format, and splitted each chromosome file to be in 2Mb segments 

5. Run the perl script '3\_Blast\_soybean\_snps\_to\_Medicagotruncatula\_genome\_to\_infer\_shared\_snps.pl'

		$ nohup perl 3_Blast_soybean_snps_to_Medicagotruncatula_genome_to_infer_shared_snps.pl 1 1 0 20 &