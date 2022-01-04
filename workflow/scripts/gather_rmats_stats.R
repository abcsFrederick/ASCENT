#!/usr/bin/env Rscript

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

# summary files only

glob_pattern="summary.txt"
sort_column="SignificantEventsJCEC"


files=list.files(pattern = glob_pattern,recursive = TRUE)
files
# functions
get_contrast<-function(fp){
  return(unlist(strsplit(fp,split="/"))[1])
}
read_input_file<-function(fp,contrast_name){
  df=read.csv(fp,header=TRUE,sep="\t",
           check.names = FALSE,
           comment.char = "#",
           strip.white = TRUE)
  df$contrast_name=contrast_name
  return(df)
}


contrasts=unlist(lapply(files, get_contrast))

for (i in 1:length(contrasts)){
  x=read_input_file(files[i],contrasts[i])
  if (i==1){
    d=x
  } else {
    d=rbind(d,x)
  }
  rm(x)
}
rm(i)

d <- arrange(d,desc(SignificantEventsJCEC))
write.table(d,file="summarys.tsv",sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)


# loop througth ASE events

event_types=c("RI","A3SS","A5SS","MXE","SE")

for (event_type in event_types){
glob_pattern=paste(event_type,"MATS.JCEC.txt",sep=".")
# sort_column="SignificantEventsJCEC"


files=list.files(pattern = glob_pattern,recursive = TRUE)
files
# functions
get_contrast<-function(fp){
  return(unlist(strsplit(fp,split="/"))[1])
}
read_input_file<-function(fp,contrast_name){
  print("Reading file")
  print(fp)
  df=read.csv(fp,header=TRUE,sep="\t",
           check.names = FALSE,
           comment.char = "#",
           strip.white = TRUE)
  if (nrow(df)!=0){
    df$contrast_name=contrast_name
  }
  return(df)
}



contrasts=unlist(lapply(files, get_contrast))

for (i in 1:length(contrasts)){
  x=read_input_file(files[i],contrasts[i])
  if (i==1){
    d=x
  } else {
    if (nrow(d)!=0){
      d=rbind(d,x)
    }
  }
  rm(x)
}
rm(i)

d<-d[order(d$FDR),]
head(d)
write.table(d,file=paste(event_type,"JCEC.all_contrasts.tsv",sep="."),sep = "\t",quote = FALSE,row.names = FALSE,col.names = TRUE)
rm(d)
}

read2df=function(f){
  d=read.csv(f,header=TRUE,
             sep="\t",
             check.names=FALSE,
             strip.white = TRUE,
             comment.char = "#")
  return(as.data.frame(d))
}
list_of_datasets<-list("Summary"=read2df("summarys.tsv"),
                       "RI"=read2df("RI.JCEC.all_contrasts.tsv"),
                       "SE"=read2df("SE.JCEC.all_contrasts.tsv"),
                       "A3SS"=read2df("A3SS.JCEC.all_contrasts.tsv"),
                       "A5SS"=read2df("A5SS.JCEC.all_contrasts.tsv"),
                       "MXE"=read2df("MXE.JCEC.all_contrasts.tsv"))


write.xlsx(list_of_datasets, file = "rmats.results.xlsx", overwrite = TRUE)