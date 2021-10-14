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
             files[sample]=source
     return files

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
        parentdir=join(RESULTSDIR,"{contrast}","isa","rsem"),
        samplestsv=SAMPLESTSV,
        contrast="{contrast}"
        outdir=join(RESULTSDIR,"{contrast}","isa"),
        rscript=join(SCRIPTSDIR,"isa.R")
    envmodules: TOOLS["R"]["version"]
    shell:"""
set -exuf -o pipefail
c1=$(echo {{params.contrast}}|awk -F"_vs_" '{{print $1}}')
c2=$(echo {{params.contrast}}|awk -F"_vs_" '{{print $2}}')
Rscript {params.rscript} \
 -p {params.parentdir} \
 -g {params.gtf} \
 -s {params.samplestsv} \
 -c $c1 \
 -d $c2 \
 -o {params.outdir}
"""

