# rules for isoformSwitchAnalzyeR
def get_rsem_files(wildcards):
    files=dict()
    contrast=wildcards.contrast
    x=CONTRASTSDF[CONTRASTSDF['name']==contrast]
    x=x.iloc[0]
    g1=x.group1
    g2=x.group2
    for group in [g1,g2]:
        for sample in GROUP2SAMPLES[group]:
            source=join(WORKDIR,"rsem","isoformcounts",sample,sample+".RSEM.isoform.results")
            if not os.path.exists(join(RESULTSDIR,contrast,"rsem")):
                os.mkdir(join(RESULTSDIR,contrast,"rsem"))
            if not os.path.exists(join(RESULTSDIR,contrast,"rsem",sample)):
                os.mkdir(join(RESULTSDIR,contrast,"rsem",sample))
            dest=join(RESULTSDIR,contrast,"rsem",sample,sample+".RSEM.isoform.results")
            os.link(source,dest)
            files[sample+"_original"]=source
            files[sample]=dest
    return files

def get_condition1(contrast):
    x=CONTRASTSDF[CONTRASTSDF['name']==contrast]
    x=x.iloc[0]
    g1=x.group1
    g2=x.group2
    return g1

def get_condition2(contrast):
    x=CONTRASTSDF[CONTRASTSDF['name']==contrast]
    x=x.iloc[0]
    g1=x.group1
    g2=x.group2
    return g2

rule isa:
    input:
        unpack(get_rsem_files),
    output:
        join(RESULTSDIR,"{contrast}","SplicingSummary.pdf")
    params:
        gtf=GTF,
        parentdir=join(RESULTSDIR,"{contrast}","rsem"),
        samplestsv=SAMPLESTSV,
        condition1=get_condition1("{contrast}"),
        condition2=get_condition2("{contrast}"),
        outdir=join(RESULTSDIR,"{contrast}"),
        rscript=join(SCRIPTSDIR,"isa.R")
    envmodules: TOOLS["R"]["version"]
    shell:"""
set -exuf -o pipefail
Rscript {params.rscript} \
 -p {params.parentdir} \
 -g {params.gtf} \
 -s {params.samplestsv} \
 -c {params.condition1} \
 -d {params.condition2} \
 -o {params.outdir}
"""

