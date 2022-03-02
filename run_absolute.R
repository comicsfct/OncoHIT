#!/usr/bin/env Rscript
library(ABSOLUTE)
args = commandArgs(trailingOnly=TRUE)
# CR002_M19_tumor.alignment.cns
file <- args[1]
sample <- sub("_.*","",file)

print(sample)
cnsdata <- read.table(file, sep="\t", quote = '', header = T, stringsAsFactors = F)
seg.dat.fn <- cnsdata[,c("chromosome","start","end","probes","log2")]
colnames(seg.dat.fn) <- c("Chromosome","Start","End","Num_Probes","Segment_Mean")
seg.dat.fn <- seg.dat.fn[!seg.dat.fn$Chromosome %in% c("chrX","chrY"),]
seg.dat.fn$Chromosome <- sub("chr","",seg.dat.fn$Chromosome)
write.table(seg.dat.fn, paste(sample, ".absolute.seg.tab", sep=""), sep="\t", quote = F, row.names = F)
seg.dat.fn <- paste(sample, ".absolute.seg.tab", sep="")

# parameters chosen from https://www.genepattern.org/modules/docs/ABSOLUTE
sigma.p <- 0 # starting point? number in the example...
# example is 0.02
max.sigma.h <- 0.015 # why??? example... it points to equation 6 in paper but no real idea why this value was chosen...
min.ploidy <- 0.95 # why?... assuming it needs at least 1 of each, except Y chromosome...
max.ploidy <- 10 # arbitrary number used in the example...

primary.disease <- "Colon Cancer"
platform <- "Illumina_WES"
sample.name <- sample
results.dir <- paste(sample, "_absolute_1", sep="")

max.as.seg.count <- 1500 # be careful all samples have less than this number, otherwise increase a bit
# example is 0
#max.non.clonal <- 0.05
#max.non.clonal <- 0.15
max.non.clonal <- 1
# example is 0
max.neg.genome <- 0.005

copy_num_type <- "total"

var_data <- read.table("somatic_mutations_filtered_imm_cohort.tab", sep="\t", header=T, stringsAsFactors = F)
var_data <- var_data[var_data$sample==sample,]

var_data <- var_data[, c("sample","ID","Chr","Start","Gene.refGene", "TUMOR")]
var_data <- var_data[! var_data$Chr %in% c("chrX","chrY"), ]
var_data$Chr <- sub("chr","",var_data$Chr) 
var_data$ID <- sub("rs.*","bydbSNP", var_data$ID)
var_data$ID <- sub("\\.","", var_data$ID)
  
var_data$t_ref_count <- sub(",.*","",sapply(strsplit(x = var_data$TUMOR, ':'), function (list){ list[2]}))
var_data$t_alt_count <- sub(".*,","",sapply(strsplit(x = var_data$TUMOR, ':'), function (list){ list[2]}))
var_data$TUMOR <- NULL
colnames(var_data) <- c("Tumor_Sample_Barcode","dbSNP_Val_Status","Chromosome","Start_position","Hugo_Symbol","t_ref_count","t_alt_count")
write.table(var_data, paste(sample, ".absolute.mut.tab", sep=""), sep="\t", quote = F, row.names = F)

maf.fn <- paste(sample, ".absolute.mut.tab", sep="")

results.dir <- paste(sample, "_absolute_maf_1", sep="")
min.mut.af=0.05
output.fn.base=sample

RunAbsolute(seg.dat.fn = seg.dat.fn, sigma.p = sigma.p, max.sigma.h = max.sigma.h, min.ploidy = min.ploidy, max.ploidy = max.ploidy, 
            primary.disease = primary.disease, platform = platform, sample.name = sample.name, results.dir = results.dir, 
            max.as.seg.count = max.as.seg.count, max.non.clonal = max.non.clonal, 
            max.neg.genome = max.neg.genome, copy_num_type = copy_num_type, 
            maf.fn = maf.fn, min.mut.af = min.mut.af, output.fn.base = output.fn.base, 
            verbose=FALSE)

