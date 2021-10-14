# rules for isoformSwitchAnalzyeR
def get_rsem_files(wildcards):
    print(wildcards)
    files=dict()
    contrast=wildcards.contrast
    x=CONTRASTSDF[CONTRASTSDF['name']==contrast]
    x=x.iloc[0]
    g1=x.group1
    g2=x.group2
    for group in [g1,g2]:
        for sample in GROUP2SAMPLES[group]:
            source=join(WORKDIR,"rsem","isoformcounts",sample,sample+".RSEM.isoform.results")
            files[sample]=source
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

rule isa_init:
    input:
        unpack(get_rsem_files)
    output:
        dummy=join(RESULTSDIR,"{contrast}","isa","dummy")
    params:
        parentdir=join(RESULTSDIR,"{contrast}","isa","rsem")
    shell:"""
set -exuf -o pipefail
if [ ! -d {params.parentdir} ];then
    mkdir -p {params.parentdir}
fi
for i in {input}; do
    bn=$(basename $i)
    sn=$(echo $bn|awk -F".RSEM" '{{print $1}}')
    mkdir -p {params.parentdir}/${{sn}}
    ln $i {params.parentdir}/${{sn}}/${{bn}}
done
touch {output.dummy}
"""


rule isa:
    input:
        rules.isa_init.output.dummy,
    output:
        join(RESULTSDIR,"{contrast}","isa","SplicingSummary.pdf")
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

