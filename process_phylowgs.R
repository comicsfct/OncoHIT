#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=6) {
  stop("args: input_vcf input_cnv input_absolute_segtab purity output_vcf output_cnv .n", call.=FALSE)
} else {
  input_vcf = args[1]
  input_cnv = args[2]
  input_absolute_segtab = args[3]
  purity = args[4]
  output_vcf = args[5]
  output_cnv = args[6]

}


cnv_data <- read.table(input_cnv, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
cnv_data <- cnv_data[,c("chromosome", "start", "end", "cn", "cn1", "cn2")]
cnv_data$chromosome <- sub("chr","",cnv_data$chromosome)
# Ignore X and Y chromosomes
cnv_data <- cnv_data[!cnv_data$chromosome %in% c("X","Y"),]

absolute_segtab <- read.table(input_absolute_segtab, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
absolute_segtab <- absolute_segtab[,c("Chromosome", "Start.bp", "End.bp", "cancer_cell_frac")]
colnames(absolute_segtab) <- c("chr", "start", "end", "ccf")
absolute_segtab$chr <- sub("chr","",absolute_segtab$chr)
absolute_segtab$ccf[is.na(absolute_segtab$ccf)] <- 1
absolute_segtab$cellular_prevalence <- purity * absolute_segtab$ccf;
absolute_segtab$ccf <- NULL

absolute_segtab$major_cn <- 1
absolute_segtab$minor_cn <- 1

for(segtab_row in 1:nrow(absolute_segtab)){
  for(cna_row in 1:nrow(cnv_data)){
    if( (absolute_segtab$chr[segtab_row] == cnv_data$chromosome[cna_row]) & 
        (absolute_segtab$end[segtab_row] > cnv_data$start[cna_row]) & 
        (absolute_segtab$end[segtab_row] < cnv_data$end[cna_row]) & 
        (absolute_segtab$start[segtab_row] < cnv_data$end[cna_row]) & 
        (absolute_segtab$start[segtab_row] > cnv_data$start[cna_row])){
      absolute_segtab[segtab_row, "major_cn"] <- cnv_data[cna_row,"cn1"]
      absolute_segtab[segtab_row, "minor_cn"] <- cnv_data[cna_row,"cn2"]
    }  
  }  
}

write.table(absolute_segtab, output_cnv, sep="\t", row.names=F, quote=F)


vcf_data <- read.table(input_vcf, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
colnames(vcf_data)[1] <- "CHROM"
vcf_data <- vcf_data[grep("clustered_events", vcf_data$FILTER, invert=T),]
vcf_data <- vcf_data[grep("t_lod_fstar", vcf_data$FILTER, invert=T),]
vcf_data$CHROM <- sub("#","", vcf_data$CHROM)
vcf_data <- vcf_data[!vcf_data$CHROM %in% c("X","Y"),]
vcf_data <- vcf_data[(nchar(vcf_data$ref)==1) & (nchar(vcf_data$alt)==1),]

vcf_data$refc <- sapply(
                               sapply(sapply(strsplit(vcf_data$TUMOR,":"),
                                             function(c){c[2]}),
                                      function(c){strsplit(c,",")}), 
                               function(vec){as.integer(vec[1])})
							   
vcf_data$altc <- sapply(
                               sapply(sapply(strsplit(vcf_data$TUMOR,":"),
                                             function(c){c[2]}),
                                      function(c){strsplit(c,",")}), 
                               function(vec){as.integer(vec[2])})
							   
vcf_data$vaf <- altc / (refc + altc)	

vcf_data <- vcf_data[vcf_data$vaf>0.05,] 	

vcf_data$INFO <- paste('VAF=',vcf_data$vaf.';t_alt_count=',t_alt_count$altc,';t_ref_count='.vcf_data$refc, sep="")
				   						   
write.table(vcf_data, output_cnv, sep="\t", row.names=F, quote=F)	

						   