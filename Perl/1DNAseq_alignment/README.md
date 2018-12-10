##This folder describe the process for maize DNAseq short reads alignment

1. Build BWA index (Linux command line)

		$cd ~/XXX/Refs/AGPV4
		$wget ftp://ftp.ensemblgenomes.org/pub/plants/release-41/fasta/zea_mays/dna/Zea_mays.AGPv4.dna.toplevel.fa.gz &
		$gunzip Zea_mays.AGPv4.dna.toplevel.fa.gz
		$ln -s Zea_mays.AGPv4.dna.toplevel.fa ZmB73_V4_bwa
		$bwa index -a bwtsw ZmB73_V4_bwa

2. Build index for fasta file

		$mv Zea_mays.AGPv4.dna.toplevel.fa ZmB73_V4_all.fa
		$samtools faidx ZmB73_V4_all.fa

3. Create sequence directory

		$samtools dict ZmB73_V4_all.fa -o ZmB73_V4_all_sam.dict


	Extract the sam header file

		$cut -f-3 ZmB73_V4_all_sam.dict > ZmB73_V4_all_sam_header

4. Prepare the SRR list, example 'sraRun\_by\_SRRs\_list' looks like following：

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

5. Run the perl script named '1\_BWA\_MEM\_run\_maize\_ATR.pl'

		$perl 1_BWA_MEM_run_maize_ATR.pl 0 4 &