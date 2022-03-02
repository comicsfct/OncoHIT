
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=4) {
  stop("args:  mutect_somatic_vcf cnvkit_cns loci output .n", call.=FALSE)
} else {
  vcf = args[1]
  cnv = args[2]
  loci = args[3]
  output = args[4]
}


vcf_data <- read.table(vcf, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
colnames(vcf_data)[1] <- "CHROM"

# only include those with PASS ? maybe also include some with traces of germline?
#vcf_data <- vcf_data[vcf_data$FILTER=="PASS",]
vcf_data <- vcf_data[grep("clustered_events", vcf_data$FILTER, invert=T),]
vcf_data <- vcf_data[grep("t_lod_fstar", vcf_data$FILTER, invert=T),]

# do not include sex chromosomes just to be safe...
vcf_data <- vcf_data[! vcf_data$CHROM %in% c("chrX", "chrY", "X", "Y"),]



cnv_data <- read.table(cnv, sep="\t", header = T)
cnv_data <- cnv_data[! cnv_data$chromosome %in% c("chrX", "chrY", "X", "Y") , ]

cnv_affected <- cnv_data[ (cnv_data$cn1 != 1) | (cnv_data$cn2 != 1),]
# more conservative... 
cnv_affected2 <- cnv_data[ (cnv_data$cn != 2),]



loci_data <-  read.table(loci, sep="\t", header=F)
colnames(loci_data)[1:3] <- c("chr","start","end")

loci_data <- loci_data[! (loci_data$chr %in% c("X","Y","chrX", "chrY")),]
# In MB
total_len <- sum(loci_data$end - loci_data$start) / 1000000

# TMB: number of mutations / MB of covered region 
tmb <- nrow(vcf_data) / total_len

# CNB: Percentage of the measured genome affected by copy number alterations...
cnb <- sum(cnv_affected$end - cnv_affected$start, na.rm = T) / sum(cnv_data$end - cnv_data$start, na.rm = T)

# CNB: Percentage of the measured genome affected by copy number alterations...
cnb2 <- sum(cnv_affected2$end - cnv_affected2$start, na.rm = T) / sum(cnv_data$end - cnv_data$start, na.rm = T)


# Return these values...

outdata <- data.frame("TMB"=signif(tmb,3),"CNB"=signif(cnb,3), "CNB_strict"=signif(cnb2,3))
write.table(outdata, output, quote=F, row.names=F)
