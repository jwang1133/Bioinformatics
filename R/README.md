## This folder contains some R scripts to show my coding experience in R

**1. Rscript '1Simple\_linear\_model\_GWAS' is for linear model of Genome-Wide Association Study on a large sample size (~1800) and large genotpe data (~21 million SNPs).**
	

- In this script, we read into a small segment from large genotype file because of the limited resource
- It facilitate parallel computing by passing parameters of chromosome through command line
- Input file format:

	**a,genotype data is in .hmp format** 

	**b,phenotype file**
 
			Taxa	PH
			L1	156.5
			L2	176
			L3	169
			L4	145
			L5	132
			L6	129
			
	**c, snpnum file**

			CHROM   POS
			1       10045
			1       10097
			1       10638
			1       10760
			1       11250
			1       11366
			1       11399
			1       11414
			1       11454

- Run script using the following command line 
		
		Rscript 1Simple_linear_model_GWAS.r 1 10 &


**2. Rscript '2Process\_multiple\_growth\_stage\_data\_and\_visualize\_data' is for processing data from multi growth stage measurement and visualizing the data**

- This script first process and extract data from two different datasets, then merge the data, further manipulate the data, and then visualize the data
- Two input files looks like:

	 **a, data1**
		
		Block	Gen	Year	NDVI_06132017	NDVI_06202017	NDVI_07062017	NDVI_07192017	NDVI_08302017
		1	1	2017	123.58	128.77	127.55	127.22	127.54
		1	119	2017	117.08	123.7	124.87	131.82	130.46
		1	232	2017	122.65	126.81	128.6	127.08	121.76
		1	63	2017	119.36	124.14	126.55	129.09	126.33
		1	29	2017	115.87	124.71	123.85	123.2	130.87
		1	65	2017	117.42	123.56	125.45	123.66	114.59
		1	111	2017	120.13	126.22	126.09	124.94	128.09
		1	17	2017	121.33	126.22	126.43	122.89	122.45
		1	18	2017	120.11	125.66	126.74	122.9	120.05
		1	139	2017	112.01	122.27	126.19	124.98	129.46

	**b, data2**
		
		Tier	Range	Pass	Block	Gen	Year	Trait	Value
		2	20	2	1	2	2015	1	70
		2	20	3	1	3	2015	1	65
		2	20	4	1	4	2015	1	80
		2	20	12	1	8	2015	1	83
		2	20	17	1	9	2017	1	.
		2	20	18	1	10	2017	2	105
		2	20	19	1	11	2017	2	134
		2	20	20	1	12	2017	2	102
		2	20	25	1	13	2017	2	220
		2	20	27	1	15	2017	2	.


**3. Rscript '3Statistical\_modeling\_of\_chromosome\_size.r' is for modeling the standardized chromosome size and the chromosome index in eukaryotic genomes with cubic distribution and gamma distribution, and plot the results out**

- In this script, we try to modeling the standardized chromosome size and the chromosome index with cubic function and gamma distribution. Parameters are estimated with non-linear least square
- Input file format (here shows the example of one genome assembly, the original data contains 913 genomes):

		Assacc	Group	SubGroup	Organism	Chr_number	Chr_length	Total_chr_nums	Genome_length	Average_chr_length	Chr_index	Standard_chr_len
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	I	201975	16	12413080	775817.5	0.031	0.26
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	VI	243427	16	12413080	775817.5	0.094	0.314
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	III	312278	16	12413080	775817.5	0.156	0.403
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	IX	422409	16	12413080	775817.5	0.219	0.544
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	VIII	531223	16	12413080	775817.5	0.281	0.685
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	V	556740	16	12413080	775817.5	0.344	0.718
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	XI	662650	16	12413080	775817.5	0.406	0.854
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	X	712875	16	12413080	775817.5	0.469	0.919
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	XIV	756632	16	12413080	775817.5	0.531	0.975
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	II	793078	16	12413080	775817.5	0.594	1.022
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	XIII	897446	16	12413080	775817.5	0.656	1.157
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	XVI	912629	16	12413080	775817.5	0.719	1.176
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	VII	1048774	16	12413080	775817.5	0.781	1.352
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	XV	1056367	16	12413080	775817.5	0.844	1.362
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	IV	1468382	16	12413080	775817.5	0.906	1.893
		GCA_000976785.2	Fungi	Ascomycetes	Saccharomyces cerevisiae YJM1202	XII	1836195	16	12413080	775817.5	0.969	2.367
		
- Fig2 illustrate how the output figure looks like
