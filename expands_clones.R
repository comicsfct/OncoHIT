library(expands)

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("args:  mutect_somatic_vcf cnvkit_cns output .n", call.=FALSE)
} else {
  vcf = args[1]
  cnv = args[2]
  output = args[3]
}

#print(vcf)
#print(cnv)

# vcf cannot have headers with ##
# Maybe better to preprocess to be sure...
vcf_data <- read.table(vcf, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
colnames(vcf_data)[1] <- "CHROM"

# only include those with PASS ? maybe also incluse some with traces of germline?
#vcf_data <- vcf_data[vcf_data$FILTER=="PASS",]
vcf_data <- vcf_data[grep("clustered_events", vcf_data$FILTER, invert=T),]
vcf_data <- vcf_data[grep("t_lod_fstar", vcf_data$FILTER, invert=T),]

# do not include sex chromosomes just to be safe...
vcf_data <- vcf_data[! vcf_data$CHROM %in% c("chrX", "chrY", "X", "Y"),]

dmcase_snv <- data.frame(chr=vcf_data$CHROM,
                         startpos=vcf_data$POS,
                         endpos=vcf_data$POS+nchar(vcf_data$REF),
                         AF_Tumor=sapply(
                           sapply(
                             sapply(strsplit(vcf_data$TUMOR,":"),
                                    function(c){c[2]}),
                             function(c){strsplit(c,",")}), 
                           function(vec){as.integer(vec[2])/(as.integer(vec[1])+as.integer(vec[2]))}),
                         PN_B=sapply(vcf_data$FILTER, FUN = function(filter){if(filter=="PASS"){return(0)} else {return(1)} }))
dmcase_snv$chr <- as.integer(sub("chr","",dmcase_snv$chr))

cnvkit_data <- read.table(cnv, sep="\t", header = T)
dmcase_cnv <- cnvkit_data[! cnvkit_data$chromosome %in% c("chrX", "chrY", "X", "Y") ,c("chromosome","start","end","cn")]
colnames(dmcase_cnv) <- c("chr","startpos","endpos", "CN_Estimate")
dmcase_cnv$chr <- as.integer(sub("chr","",dmcase_cnv$chr))

dm <- assignQuantityToMutation(as.matrix(dmcase_snv),as.matrix(dmcase_cnv),"CN_Estimate")

# Eventually take these parameters as input...
max_PM=6; maxS=0.7; precision=0.018; plotF=1;
set.seed(101) # Otherwise it won't be deterministic...
cfd <- computeCellFrequencyDistributions(dm, max_PM=max_PM, p=precision, nc=1)
toUseIdx <- which(apply(is.finite(cfd$densities),1,all) )
SPs <- clusterCellFrequencies(cfd$densities[toUseIdx,], p=precision)

SPs <- SPs[SPs[,"score"] <= maxS,]

print(SPs)
#add line: if file missing
write.table(SPs,output,sep="\t",quote=F)


