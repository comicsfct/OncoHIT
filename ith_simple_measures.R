
args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("args:  mutect_somatic_vcf cnvkit_cns output .n", call.=FALSE)
} else {
  vcf = args[1]
  cnv = args[2]
  output = args[3]
}

vcf_data <- read.table(vcf, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
colnames(vcf_data)[1] <- "CHROM"

# only include those with PASS? (this is likely too strict) 
# maybe also include some with traces of germline?
#vcf_data <- vcf_data[vcf_data$FILTER=="PASS",]
vcf_data <- vcf_data[grep("clustered_events", vcf_data$FILTER, invert=T),]
vcf_data <- vcf_data[grep("t_lod_fstar", vcf_data$FILTER, invert=T),]

# do not include sex chromosomes just to be safe...
vcf_data <- vcf_data[! vcf_data$CHROM %in% c("chrX", "chrY", "X", "Y"),]

dmcase_snv <- data.frame(chromosome=vcf_data$CHROM,
                         start=vcf_data$POS,
                         end=vcf_data$POS+nchar(vcf_data$REF),
                         AF_Tumor=sapply(
                           sapply(
                             sapply(strsplit(vcf_data$TUMOR,":"),
                                    function(c){c[2]}),
                             function(c){strsplit(c,",")}), 
                           function(vec){as.integer(vec[2])/(as.integer(vec[1])+as.integer(vec[2]))}),
						   stringsAsFactors=F)

cnvkit_data <- read.table(cnv, sep="\t", header = T,stringsAsFactors=F)
dmcase_cnv <- cnvkit_data[! (cnvkit_data$chromosome %in% c("chrX", "chrY", "X", "Y")),]

dmcase_cnv <- dmcase_cnv[(dmcase_cnv$cn1 != 1) | (dmcase_cnv$cn2 != 1),]
dmcase_cnv <- dmcase_cnv[!is.na(dmcase_cnv$chromosome),]
# More conservative
#dmcase_cnv <- dmcase_cnv[dmcase_cnv$cn!=2,]

# remove those that intersect
# very inefficient, but should be ok as there should be relatively few (somatic) SNVs and CNAs
remove <- c()
for(snp_row in 1:nrow(dmcase_snv)){
  for(cna_row in 1:nrow(dmcase_cnv)){
    if( (dmcase_snv$chromosome[snp_row] == dmcase_cnv$chromosome[cna_row]) & 
        (dmcase_snv$end[snp_row] > dmcase_cnv$start[cna_row]) & 
        (dmcase_snv$start[snp_row] < dmcase_cnv$end[cna_row])){
      remove <- c(snp_row,remove)
    }  
  }  
}
# Either use an apply, or use a specialized package
#apply( X = dmcase_snv_nocna, MARGIN = 1, FUN = function(row){} )
dmcase_snv_nocna <- dmcase_snv
if(length(remove)>0){
  dmcase_snv_nocna <- dmcase_snv[-remove,]
}

# Give the two numbers, with all and without those potentially affected by CNAs

# MATH: 
math <- 100 * mad(dmcase_snv$AF_Tumor) / median(dmcase_snv$AF_Tumor)
math_nocna <- 100 * mad(dmcase_snv_nocna$AF_Tumor) / median(dmcase_snv_nocna$AF_Tumor)

# TH Index
dmcase_snv$bin <- floor(10*dmcase_snv$AF_Tumor)
dmcase_snv_nocna$bin <- floor(10*dmcase_snv_nocna$AF_Tumor)

th <- sum(sapply(table(dmcase_snv$bin) / nrow(dmcase_snv), FUN = function(p){ -1 * p*log(p)} ))
th_nocna <- sum(sapply(table(dmcase_snv_nocna$bin) / nrow(dmcase_snv_nocna), FUN = function(p){ -1 * p*log(p)} ))

# Return these values...

outdata <- data.frame("MATH"=signif(math,3),"MATH_noCNA"=signif(math_nocna,3),"TH"=signif(th,3),"TH_noCNA"=signif(th_nocna,3))
write.table(outdata, output, quote=F, row.names=F)
