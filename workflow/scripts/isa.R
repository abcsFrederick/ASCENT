#!/usr/bin/env Rscript

.libPaths( c( "/data/Ziegelbauer_lab/rlib",.libPaths() ) )
suppressPackageStartupMessages(library("argparse"))

# create parser object
parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 

parser$add_argument("-p", "--parentdir", 
                    type="character", 
                    help="parent directory of RSEM isoform counts (1 folder per sample)",
                    required=TRUE)
parser$add_argument("-g", "--gtf", 
                    type="character", 
                    help="path to the GTF annotations file",
                    required=TRUE)
parser$add_argument("-s", "--samplestsv", 
                    type="character",
                    help="samples.tsv file with group column which will be used as condition",
                    required=TRUE)
parser$add_argument("-c", "--condition1", 
                    type="character",
                    help="comma separated list of condition 1's.... condition1 vs condition2 contrast will be made for each pair",
                    required=TRUE)
parser$add_argument("-d", "--condition2", 
                    type="character",
                    help="comma separated list of condition 2's.... condition1 vs condition2 contrast will be made for each pair",
                    required=TRUE)
parser$add_argument("-o", "--outdir", 
                    type="character",
                    help="output dir for all output files",
                    required=TRUE)

args <- parser$parse_args()


suppressPackageStartupMessages(library("IsoformSwitchAnalyzeR"))
suppressPackageStartupMessages(library("tidyverse"))
debug=0

parentdir=args$parentdir
samplestsv=args$samplestsv
gtf=args$gtf
c1s=unlist(strsplit(args$condition1,","))
c2s=unlist(strsplit(args$condition2,","))
outdir=args$outdir

if(debug==1){
  # parentdir="/Volumes/Ziegelbauer_lab/rmats_testing2/rsem/isoformcounts"
  parentdir="/Volumes/Ziegelbauer_lab/rmats_testing2/results/isoformSwitchAnalyzer/SD1_ORF57KO-24h_vs_SD1_ORF57KO-0h/rsem"
  samplestsv="/Volumes/Ziegelbauer_lab/rmats_testing2/samples.tsv"
  gtf="/Volumes/Ziegelbauer_lab/resources/hg38/hg38.v36.gtf"
  rdata=paste0(parentdir,"/isoformSwitchAnalyzer.Rdata")
  c1s=unlist(strsplit(c("SD1_ORF57KO-24h","SD1_ORF57KO-0h"),","))
  c2s=unlist(strsplit(c("SD1_ORF57KO-0h","SD1_ORF57KO-24h"),","))
  outdir="/Volumes/Ziegelbauer_lab/rmats_testing2/results/isoformSwitchAnalyzer/SD1_ORF57KO-24h_vs_SD1_ORF57KO-0h"
}

stopifnot("Number of condition1s and conditions2s do not match!" = length(c1s)==length(c2s))

comp2make=data.frame(condition_1=c1s,condition_2=c2s)

# import RSEM results
rsemcounts=importIsoformExpression(parentDir = parentdir,
                        )
### Create switchAnalyzeRlist
samples=read.table(samplestsv,header=TRUE,sep="\t")

samples[,c("sampleName","group")] %>% 
  as.data.frame() %>%
  column_to_rownames(var = "sampleName") -> sample2group

sampleID=colnames(rsemcounts$abundance)[-1]
condition=sample2group[sampleID,]
myDesign=data.frame(sampleID=sampleID,condition=condition)


aSwitchList <- importRdata(
  isoformCountMatrix = rsemcounts$counts,
  isoformRepExpression = rsemcounts$abundance,
  isoformExonAnnoation = gtf,
  designMatrix = myDesign,
  showProgress = TRUE,
  comparisonsToMake = comp2make
  )


# filter the switchAnalyzeRlist
filteredSwitchList <- preFilter(
  switchAnalyzeRlist = aSwitchList,
  geneExpressionCutoff = 1,
  isoformExpressionCutoff = 0,
  IFcutoff=0.01,
  alpha=0.05,
  dIFcutoff = 0.1,
  quiet=FALSE
  )

# run DEXseq
# Perform test
aSwitchListAnalyzed <- isoformSwitchTestDEXSeq(
  switchAnalyzeRlist = filteredSwitchList,
  reduceToSwitchingGenes=TRUE,
  quiet = FALSE
)
# extractSwitchSummary(aSwitchListAnalyzed,onlySigIsoforms=TRUE)


merge(myDesign,aSwitchListAnalyzed$designMatrix,by=c("sampleID")) %>% 
  select(-c("sampleID")) %>%
  unique() %>%
  remove_rownames() -> conditions_old_new

conditions_old_new %>%
  column_to_rownames(var = "condition.x") -> conditions_old2new
conditions_old_new %>%
  column_to_rownames(var = "condition.y") -> conditions_new2old

# analyze alternate splicing
aSwitchListAnalyzed <- analyzeAlternativeSplicing(
  switchAnalyzeRlist = aSwitchListAnalyzed,
  quiet=FALSE
)

pdf(paste0(outdir,"/SplicingSummary.pdf"))
extractSplicingSummary(
  aSwitchListAnalyzed,
  asFractionTotal = FALSE,
  plotGenes=FALSE
)
dev.off()

pdf(paste0(outdir,"/SplicingEnrichment.pdf"))
splicingEnrichment <- extractSplicingEnrichment(
  aSwitchListAnalyzed,
  splicingToAnalyze='all',
  returnResult=TRUE,
  returnSummary=TRUE
)
dev.off()
splicingEnrichment$condition_1=conditions_new2old[splicingEnrichment$condition_1,]
splicingEnrichment$condition_2=conditions_new2old[splicingEnrichment$condition_2,]
splicingEnrichment$Comparison=NULL

pdf(paste0(outdir,"/GenomeWideSplicing.pdf"))
genomeWideSplicing<-extractSplicingGenomeWide(
  aSwitchListAnalyzed,
  featureToExtract = 'all', 
  splicingToAnalyze = 'all',
  plot=TRUE,
  returnResult=TRUE  
)
dev.off()
genomeWideSplicing %>% separate(col = "comparison", into=c("condition_1","condition_2"),sep=" vs ") -> genomeWideSplicing
genomeWideSplicing$condition_1=conditions_new2old[genomeWideSplicing$condition_1,]
genomeWideSplicing$condition_2=conditions_new2old[genomeWideSplicing$condition_2,]

# results
tsg=extractTopSwitches(
  switchAnalyzeRlist=aSwitchListAnalyzed,
  extractGenes=TRUE,
  alpha=0.05,
  dIFcutoff = 0.1,
  sortByQvals=TRUE,
  n=100000
)

tsi=extractTopSwitches(
  switchAnalyzeRlist=aSwitchListAnalyzed,
  extractGenes=FALSE,
  alpha=0.05,
  dIFcutoff = 0.1,
  sortByQvals=TRUE,
  n=100000
)
tsg$condition_1=conditions_new2old[tsg$condition_1,]
tsg$condition_2=conditions_new2old[tsg$condition_2,]
tsi$condition_1=conditions_new2old[tsi$condition_1,]
tsi$condition_2=conditions_new2old[tsi$condition_2,]
tsg$Rank=NULL
tsi$Rank=NULL

subset_and_write <- function(dt,c1,c2,of){
  x=dt[dt$condition_1==c1 & dt$condition_2==c2,]
  write.table(x,file=of,sep = "\t",quote = FALSE,row.names = FALSE)
}

for (row in 1:nrow(comp2make)){
  c1=comp2make[row,c("condition_1")]
  c2=comp2make[row,c("condition_2")]
  prefix=paste0(c1,"_vs_",c2)
  
  # tsg
  outtsv=paste0(outdir,"/",prefix,".top_switched_genes.tsv")
  subset_and_write(dt=tsg,c1=c1,c2=c2,of=outtsv)

  # tsi
  outtsv=paste0(outdir,"/",prefix,".top_switched_isoforms.tsv")
  subset_and_write(dt=tsi,c1=c1,c2=c2,of=outtsv)

  # splicingEnrichment
  outtsv=paste0(outdir,"/",prefix,".splicingEnrichment.tsv")
  subset_and_write(dt=splicingEnrichment,c1=c1,c2=c2,of=outtsv)

  # genomeWideSplicing
  outtsv=paste0(outdir,"/",prefix,".genomeWideSplicing.tsv")
  subset_and_write(dt=genomeWideSplicing,c1=c1,c2=c2,of=outtsv)
  
}




# switchPlotTopSwitches( aSwitchListAnalyzed )
# switchPlotIsoExp(aSwitchListAnalyzed, gene = 'CLPTM1L')
# switchPlotTranscript(aSwitchListAnalyzed, gene = 'CLPTM1L')
# switchPlotIsoUsage(aSwitchListAnalyzed, gene = 'CLPTM1L')

# save image for future plots etc.
save.image(file=paste0(outdir,"/isoformSwitchAnalyzer.Rdata"))
