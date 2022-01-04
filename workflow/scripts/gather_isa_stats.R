#!/usr/bin/env Rscript
# This script is used to concatenate (row bind) all results from IsoformSwitchAnalyzer.
# ISA generated 4 types of files:
# 1. top_switched_genes.tsv which are sorted by gene_switch_q_value
# 2. top_switched_isoforms.tsv which are sorted by isoform_switch_q_value
# 3. splicingEnrichment.tsv which are sorted by propUpQval
# 4. genomeWideSplicing.tsv which are sorted by wilcoxQval
#
# Eg:
#
# Rscript gather_isa_stats.R \
#  -d /data/Ziegelbauer_lab/ccbr1140/samples_1/results \
#  -p top_switched_genes.tsv \
#  -s gene_switch_q_value \
#  -o top_switched_genes.all_samples.tsv

#
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-d", "--dirname", 
                    type="character", 
                    help="results dir of the pipeline or dir where to search recursively for files with given pattern",
                    required=TRUE)


args <- parser$parse_args()

suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("openxlsx"))

wd=args$dirname
setwd(wd)

globpatterns=c("top_switched_genes.tsv","top_switched_isoforms.tsv","splicingEnrichment.tsv","genomeWideSplicing.tsv")
sortcolumns=c("gene_switch_q_value","isoform_switch_q_value","propUpQval","wilcoxQval")
outfiles=c("top_switched_genes.all_samples.tsv","top_switched_isoforms.all_samples.tsv","splicingEnrichment.all_samples.tsv","genomeWideSplicing.all_samples.tsv")

for (i in 1:length(globpatterns)){
  glob_pattern=globpatterns[i]
  sort_column=sortcolumns[i]
  outfile=outfiles[i]

  files=list.files(pattern = glob_pattern,recursive = TRUE)
  files

  read_input_file<-function(fp){
    df=read.csv(fp,header=TRUE,sep="\t",
                check.names = FALSE,
                comment.char = "#",
                strip.white = TRUE)
    return(df)
  }

  for (i in 1:length(files)){
    x=read_input_file(files[i])
    if (i==1){
      d=x
    } else {
      d=rbind(d,x)
    }
    rm(x)
  }
  rm(i)

  d=d[order(d[,c(sort_column)]),]
  write.table(d,
              file=outfile,
              row.names = FALSE,
              col.names = TRUE,
              quote = FALSE,
              sep="\t")
}

read2df=function(f){
  d=read.csv(f,header=TRUE,
             sep="\t",
             check.names=FALSE,
             strip.white = TRUE,
             comment.char = "#")
  return(as.data.frame(d))
}

list_of_datasets<-list("top_switched_genes"=read2df("top_switched_genes.all_samples.tsv"),
                       "top_switched_isoforms"=read2df("top_switched_isoforms.all_samples.tsv"),
                       "splicingEnrichment"=read2df("splicingEnrichment.all_samples.tsv"),
                       "genomeWideSplicing"=read2df("genomeWideSplicing.all_samples.tsv"))


write.xlsx(list_of_datasets, file = "isa.results.xlsx", overwrite = TRUE)