library(sleuth)
library(data.table)
#extract results and see most significants results
library(dplyr)

#read in the table you made describing samples and kallisto output 
stab <- read.table("input_sleuth.txt", header=TRUE, stringsAsFactors = FALSE, sep='\t')

#initialize sleuth object
so <- sleuth_prep(stab)

#differential expression analysis comparing scramble to knockdown

#fit a model comparing the two conditions
so <- sleuth_fit(so, ~condition, 'full')
#fit the reduced model to compare in the likelihood ratio test 
so <- sleuth_fit(so, ~1, 'reduced')
#perform the likelihood ratio test for differential expression between conditions 
so <- sleuth_lrt(so, 'reduced','full')


#extract the test results from the sleuth object 
sleuth_table <- sleuth_results(so, 'reduced:full','lrt',show_all=FALSE)
#filter most significant results (FDR/qval < 0.05) and sort by pval
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05) %>% dplyr::arrange(pval)

sig_sleuth <- sleuth_significant %>% select(target_id, test_stat, pval, qval)
#write target id, test stat, pval and qval for significant transcript
#include header, tab-delimit
write.table(sig_sleuth, file="R_sleuth_output.txt",quote= FALSE,row.names= FALSE,sep = "\t")