import sys
import os
import pandas as pd
import yaml
from snakemake.utils import validate
import subprocess


def get_file_size(filename):
    filename=filename.strip()
    if check_readaccess(filename):
        return os.stat(filename).st_size

def get_fastqs(wildcards):
    d=dict()
    d["R1"]=SAMPLESDF["R1"][wildcards.sample]
    d["R2"]=SAMPLESDF["R2"][wildcards.sample]
    return d

def get_peorse(wildcards):
    return SAMPLESDF["PEorSE"][wildcards.sample]

def check_existence(filename):
    """Checks if file exists on filesystem
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    if not os.path.exists(filename):
        sys.exit("File: {} does not exists!".format(filename))
    return True


def check_readaccess(filename):
    """Checks permissions to see if user can read a file
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    check_existence(filename)
    if not os.access(filename,os.R_OK):
        sys.exit("File: {} exists, but user cannot read from file due to permissions!".format(filename))
    return True


def check_writeaccess(filename):
    """Checks permissions to see if user can write to a file
    :param filename <str>: Name of file to check
    """
    filename=filename.strip()
    check_existence(filename)
    if not os.access(filename,os.W_OK):
        sys.exit("File: {} exists, but user cannot write to file due to permissions!".format(filename))
    return True

##### load config and sample sheets #####

#validate(config, "config.schema.yaml")

## set memory limit 
## used for sambamba sort, etc
MEMORYG="100G"

#resouce absolute path
WORKDIR=config['workdir']
CONFIGFILE=join(WORKDIR,"config.yaml")
check_readaccess(CONFIGFILE)
configfile: CONFIGFILE
with open(CONFIGFILE) as f:
    CONFIG = yaml.safe_load(f)

SCRIPTSDIR=config['scriptsdir']
RESOURCESDIR=config['resourcesdir']
if not os.path.exists(join(WORKDIR,"fastqs")):
    os.mkdir(join(WORKDIR,"fastqs"))

for f in ["samples", "tools", "cluster"]:
    check_readaccess(config[f])

SAMPLESDF = pd.read_csv(config["samples"],sep="\t",header=0,index_col="sampleName")
SAMPLES = list(SAMPLESDF.index)
SAMPLESDF["R1"]=join(RESOURCESDIR,"dummy")
SAMPLESDF["R2"]=join(RESOURCESDIR,"dummy")
SAMPLESDF["PEorSE"]="PE"

GROUPS=list(set(list(SAMPLESDF["group"])))
GROUP2SAMPLES=dict()
for g in GROUPS:
    GROUP2SAMPLES[g]=list(SAMPLESDF[SAMPLESDF.group==g].index)

for sample in SAMPLES:
    R1file=SAMPLESDF["path_to_R1_fastq"][sample]
    R2file=SAMPLESDF["path_to_R2_fastq"][sample]
    # print(sample,R1file,R2file)
    check_readaccess(R1file)
    R1filenewname=join(WORKDIR,"fastqs",sample+".R1.fastq.gz")
    if not os.path.exists(R1filenewname):
        os.symlink(R1file,R1filenewname)
        # os.link(R1file,R1filenewname)
    SAMPLESDF.loc[[sample],"R1"]=R1filenewname
    if str(R2file)!='nan':
        check_readaccess(R2file)
        R2filenewname=join(WORKDIR,"fastqs",sample+".R2.fastq.gz")
        if not os.path.exists(R2filenewname):
            os.symlink(R2file,R2filenewname)
            # os.link(R2file,R2filenewname)
        SAMPLESDF.loc[[sample],"R2"]=R2filenewname
    else:
        SAMPLESDF.loc[[sample],"PEorSE"]="SE"
        sys.exit(" This pipeline only works for PAIRED end samples. If your samples are paired end, then please verify that they have been entered correctly in samples.tsv ")
# print(SAMPLESDF)
# sys.exit()

#########################################################
# READ IN TOOLS REQUIRED BY PIPELINE
# THESE INCLUDE LIST OF BIOWULF MODULES (AND THEIR VERSIONS)
# MAY BE EMPTY IF ALL TOOLS ARE DOCKERIZED
#########################################################
## Load tools from YAML file
try:
    TOOLSYAML = config["tools"]
except KeyError:
    TOOLSYAML = join(WORKDIR,"tools.yaml")
check_readaccess(TOOLSYAML)
with open(TOOLSYAML) as f:
    TOOLS = yaml.safe_load(f)
#########################################################



#########################################################
# READ CLUSTER PER-RULE REQUIREMENTS
#########################################################

## Load cluster.json
try:
    CLUSTERJSON = config["cluster"]
except KeyError:
    CLUSTERJSON = join(WORKDIR,"cluster.json")
check_readaccess(CLUSTERJSON)
with open(CLUSTERJSON) as json_file:
    CLUSTER = json.load(json_file)

## Create lambda functions to allow a way to insert read-in values
## as rule directives
getthreads=lambda rname:int(CLUSTER[rname]["threads"]) if rname in CLUSTER and "threads" in CLUSTER[rname] else int(CLUSTER["__default__"]["threads"])
getmemg=lambda rname:CLUSTER[rname]["mem"] if rname in CLUSTER else CLUSTER["__default__"]["mem"]
getmemG=lambda rname:getmemg(rname).replace("g","G")
#########################################################

#########################################################
# SET OTHER PIPELINE GLOBAL VARIABLES
#########################################################


RESULTSDIR=join(WORKDIR,"results")
if not os.path.exists(RESULTSDIR):
    os.mkdir(RESULTSDIR)
QCDIR=join(RESULTSDIR,"QC")

print("# Pipeline Parameters:")
print("#"*100)
print("# Working dir :",WORKDIR)
print("# Results dir :",RESULTSDIR)
print("# Scripts dir :",SCRIPTSDIR)
print("# Resources dir :",RESOURCESDIR)
print("# Cluster json :",CLUSTERJSON)


GTF=config['ref_gtf']
check_readaccess(GTF)
print("# Annotations file :",GTF)

REFFA=config['ref_fa']
check_readaccess(REFFA)
print("# Reference fasta :",REFFA)

STARINDEXDIR=join(WORKDIR,"star_index_dir")
if not os.path.exists(STARINDEXDIR):
    os.mkdir(STARINDEXDIR)
print("# STAR index dir  :",STARINDEXDIR)

#########################################################
# get maxrl
#########################################################
filelist=""
for r1 in SAMPLESDF["R1"].tolist():
    filelist=filelist+" "+r1

cmd=SCRIPTSDIR+"/_printmaxrl2 "+filelist

MAXRL=subprocess.check_output(cmd, stderr=subprocess.STDOUT, shell = True)
MAXRL=int(MAXRL.strip())

CONTRASTSDF=pd.read_csv(config['contrasts'],sep="\t",names=["group1","group2"])
CONTRASTSDF["name"] = CONTRASTSDF["group1"] + "_vs_" + CONTRASTSDF["group2"]
CONTRASTS=CONTRASTSDF["name"].to_list()
g1s=CONTRASTSDF['group1'].to_list()
g2s=CONTRASTSDF['group2'].to_list()
GROUPSINCONTRASTS=g1s+g2s
GROUPSINCONTRASTS=list(set(GROUPSINCONTRASTS))


GROUP2RMATSTXT=dict()
RMATSTXTDIR=join(WORKDIR,"rmats_txt_files")
BAMDIR=join(WORKDIR,"bam")
if not os.path.exists(BAMDIR):
    os.mkdir(BAMDIR)
if not os.path.exists(RMATSTXTDIR):
    os.mkdir(RMATSTXTDIR)
for group in GROUPSINCONTRASTS:
    group_txt=[]
    for sample in GROUP2SAMPLES[group]:
        bamfile=join(BAMDIR,sample+".Aligned.sortedByCoord.out.bam")
        group_txt.append(bamfile)
    group_txt=",".join(map(str,group_txt))
    outfilename=join(RMATSTXTDIR,group+".txt")
    GROUP2RMATSTXT[group]=outfilename
    outfile=open(outfilename,'w')
    outfile.write(group_txt)
    outfile.close()
