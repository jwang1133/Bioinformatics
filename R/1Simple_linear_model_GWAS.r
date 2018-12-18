###The purpose of this script:
###1. To perform Genome-Wide Association Analysis with simple linear model on a population with ~1800 lines, 21M SNPs and on plant height traits
###2. Because of the limited computing resources we have, we read into small segments from each chromosome to facilitate the computing, this also enable the parallel computing

##Pass the chromosome information from command line
Args <- commandArgs(TRUE)
ch_s <- as.numeric(Args[1])
ch_e <- as.numeric(Args[2])

##Read into the phenotype file, and put the phenotype into a vector
Pheno <- read.table('/XXX/pheno/AmesDP_pheno', header = T, sep = "\t")
Y_ori  <- as.vector(Pheno$PH)

##create a vector which store the snp number for each of the chromosome
snpnum <- read.table('/XXX/geno/AmesDP_snpnum', header = T, sep = "\t")
SNP_num <- rep (0, 10)
for (i in 1:10){
  data <- snpnum[snpnum$CHROM==i, ] 
  SNP_num[i] <- nrow(data)
}

##define the size of each segments
n_rows <- 10000

for (ch in ch_s:ch_e) {
 ch_hmp_file <- paste('/XXX/geno/AmesDP_ch', ch,  '.hmp',   sep = '')
##split the genotype of each chromsome into segments
 snp_num <- SNP_num[ch]
 L <- floor(snp_num/n_rows) + 1
 
 for (j in 1:L){

  n_s <- (j - 1) * n_rows  
  n_e <- j * n_rows + 1
  ch_hmp <- read.table(ch_hmp_file, header = T, sep = "\t" , skip = n_s, nrows = n_rows)
  
  output_dir <- '/XXX/'
  if (!file.exists(output_dir)) {dir.create(output_dir, recursive = T) }
  
  simple_GWAS <- paste(output_dir, group, '_simple_ch', ch,'.hmp', sep = '')
  
  N <- nrow(ch_hmp)
  M <- nrow(Pheno)
  P_values <- c()
  F_values <- c()
  
  ##create dataframe to store the GWAS results
  result_df <- ch_hmp[,c(1:4)]
  ##loop through each SNP to perform the test
  for (i in 1:N) {
  
    all_samples <- as.vector(as.matrix(ch_hmp[i,]))
    ##remove the first 11 value in the SNP vector 
    X_ori <- all_samples[-c(1:11)]
    X <- X_ori
    Y <- Y_ori
    bad_samples <- which(X_ori == "N")
    
    ##remove missing values
    if (length(bad_samples) > 0) {
      X <- X_ori[-bad_samples];
      Y <- Y_ori[-bad_samples];
    }
    
    ##Run linear model for each of the SNP
    X <- as.factor(X);
    Y_X_lm <- lm(Y ~ X);
    P <- -log10(anova(Y_X_lm)[1,5])
    F <- anova(Y_X_lm)[1,4]
    P_values[i] <- P
    F_values[i] <- F
  }
  
   result_df$P <- P_values
   result_df$F <- F_values
   
   ##write out the results
   write.table(result_df, file = simple_GWAS, append = T, sep = "\t", col.names=FALSE)
  }
}
