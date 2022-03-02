#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

if (length(args)!=3) {
  stop("args:  sample_list output type.n", call.=FALSE)
} else {
  sample_list = args[1]
  output = args[2]
  type = args[3]
}

#print(sample_list)
#print(output)
# Make it immediately generic...

# Test1 ... 
# remove the ## header sections from the VCF
#vcf = "angiogenesis/snvs/49_ctrl1.alignment.bam_VS_MDA_MB_231.alignment.bam.intervals.bed.vcf.nohead"
#cnv = "angiogenesis/cnv_kit_test2/49_ctrl1.alignment.call.allele.cns"

#sample_list <- "samples.tab"

sample_data <- read.table(sample_list, sep="\t", header=F, stringsAsFactors=F)
colnames(sample_data) <- c("sample_id","vcf","cns")

full_data <- list()
unique_mutation_ids <- c()
for(lineid in 1:nrow(sample_data)){
  sample_id <- sample_data$sample_id[lineid]
  vcf <- sample_data$vcf[lineid]
  #print(vcf)

  vcf_data <- read.table(vcf, sep="\t", header = T, quote = '', comment.char = '', stringsAsFactors=F)
  colnames(vcf_data)[1] <- "CHROM"
  
  pyclone_data <- data.frame(mutation_id=paste(vcf_data$CHROM,vcf_data$POS, sep="_"), 
                             sample_id = sample_id,
                             chr=vcf_data$CHROM,
                             start=vcf_data$POS,
                             end=vcf_data$POS+nchar(vcf_data$REF),
                             ref_counts=sapply(
                               sapply(sapply(strsplit(vcf_data$TUMOR,":"),
                                             function(c){c[2]}),
                                      function(c){strsplit(c,",")}), 
                               function(vec){as.integer(vec[1])}),
                             alt_counts=sapply(
                               sapply(sapply(strsplit(vcf_data$TUMOR,":"),
                                             function(c){c[2]}),
                                      function(c){strsplit(c,",")}), 
                               function(vec){as.integer(vec[2])}),
                             filter=vcf_data$FILTER,
                             major_cn=0,
                             minor_cn=0,
                             normal_cn=2,
			     stringsAsFactors = F)
				 
  # Remove the most likely artifacts...				 
  pyclone_data <- pyclone_data[grep("clustered_events", pyclone_data$filter, invert=T),]
  pyclone_data <- pyclone_data[grep("t_lod_fstar", pyclone_data$filter, invert=T),]

  #unique_mutation_ids <- unique(c(unique_mutation_ids, pyclone_data$mutation_id[pyclone_data$filter=="PASS"]))
  unique_mutation_ids <- unique(c(unique_mutation_ids, pyclone_data$mutation_id))
  
  full_data[[sample_id]][["snv"]] <- pyclone_data
  
  cns <- sample_data$cns[lineid]
  #print(cns) 
  full_data[[sample_id]][["cnv"]] <- read.table(cns, sep="\t", header = T)
  
}


# TYPE: NULL (alt=ref=0); AVG (ref=AVG(total),alt=0); BAM (read from bam files) 
#type <- "NULL"

full_pyclone_data <- data.frame()

for(sample in names(full_data)){
  
  sample_pyclone_data <- full_data[[sample]][["snv"]]
  # only keep those that matter... namely the ones with PASS
  sample_pyclone_data <- sample_pyclone_data[sample_pyclone_data$mutation_id %in% unique_mutation_ids,]
  
  # Now add those that are in other samples and not in this one
  positions_to_add <- unique_mutation_ids[! (unique_mutation_ids %in% sample_pyclone_data$mutation_id)]
  
  if(length(positions_to_add) > 0){
  
    sample_pyclone_data_to_add <- data.frame(mutation_id=positions_to_add, 
                                           sample_id = sample,
                                           chr=sub("_.*","",positions_to_add),
                                           start=as.numeric(sub(".*_","",positions_to_add)),
                                           # simplification
                                           end=as.numeric(sub(".*_","",positions_to_add))+1,
                                           # type=NULL
                                           ref_counts=0, 
                                           alt_counts=0,
                                           # Not actually relevant the value
                                           filter="EXTRA",
                                           major_cn=0,
                                           minor_cn=0,
                                           normal_cn=2,
					   stringsAsFactors = F)
  
    if(type=="AVG"){
      average <- mean(sample_pyclone_data$ref_counts + sample_pyclone_data$alt_counts)
      sample_pyclone_data_to_add$ref_counts <- trunc(average)
    }
  
    # if(type=="BAM"){
    #   # READ data from BAM file using RSamtools
    # }
  
    # Make sure the order is the same
    sample_pyclone_data <- rbind(sample_pyclone_data,sample_pyclone_data_to_add)
  }	
  sample_pyclone_data <- sample_pyclone_data[order(sample_pyclone_data$mutation_id),]
  
  cnv_data <- full_data[[sample]][["cnv"]]
  
  # append info
  # very inefficient, but should be ok as there should be relatively few (somatic) SNVs and CNAs
  for(snp_row in 1:nrow(sample_pyclone_data)){
    for(cna_row in 1:nrow(cnv_data)){
      #Force it to be contained and not just intersecting though it won't make much difference...
      if( (sample_pyclone_data$chr[snp_row] == cnv_data$chr[cna_row]) & 
          (sample_pyclone_data$end[snp_row] > cnv_data$start[cna_row]) & 
          (sample_pyclone_data$end[snp_row] < cnv_data$end[cna_row]) & 
          (sample_pyclone_data$start[snp_row] < cnv_data$end[cna_row]) & 
          (sample_pyclone_data$start[snp_row] > cnv_data$start[cna_row])){
        sample_pyclone_data[snp_row, "major_cn"] <- cnv_data[cna_row,"cn1"]
        sample_pyclone_data[snp_row, "minor_cn"] <- cnv_data[cna_row,"cn2"]
      }  
    }  
  }
  
  full_pyclone_data <- rbind(full_pyclone_data, sample_pyclone_data)
  
}

# Check
sum(table(full_pyclone_data$sample_id)==length(unique_mutation_ids))==length(names(full_data))


# Sometimes the CNA allelic-specific copy number comes with NA, so need to filter those out
# Pyclone then ignores mutations for which there is no data for all samples... 
full_pyclone_data <- full_pyclone_data[complete.cases(full_pyclone_data),]

# only include those with PASS ? maybe also include some with traces of germline?
# Filter was already applied when defining unique mutation_ids 
#pyclone_data <- pyclone_data[pyclone_data$filter=="PASS",]

# do not include sex chromosomes just to be safe... otherwise add as extra column: male/female
full_pyclone_data <- full_pyclone_data[! full_pyclone_data$chr %in% c("chrX", "chrY", "X", "Y"),]

full_pyclone_data <- full_pyclone_data[,c("mutation_id", "sample_id", "ref_counts", "alt_counts", "major_cn", "minor_cn", "normal_cn")]

write.table(full_pyclone_data, output, sep="\t", quote = F, row.names = F)
