rule create_rsem_index:
    input:
        reffa=REFFA,
        gtf=GTF
    output:
        ump=join(RSEMINDEXDIR,"ref.transcripts.ump")
    params:
        rsemindexdir=RSEMINDEXDIR,
    envmodules: TOOLS['rsem']['version']
    threads: getthreads("create_rsem_index")
    shell:"""
set -euf -o pipefail
cd {params.rsemindexdir}
rsem-prepare-reference -p {threads} --gtf {input.gtf} {input.reffa} ref
rsem-generate-ngvector ref.transcripts.fa ref.transcripts
"""

rule create_bed12:
    input:
        gtf=GTF
    output:
        bed12=join(RSEMINDEXDIR,"ref.bed12")
    params:
        rsemindexdir=RSEMINDEXDIR
    envmodules:
        TOOLS['ucsc']['version'], 
    shell:"""
set -euf -o pipefail
cd {params.rsemindexdir}
genesgtf={input.gtf}
bn=$(basename $genesgtf)
gtfToGenePred -genePredExt -geneNameAsName2 $genesgtf ${{bn}}.tmp
awk -v OFS="\\\t" '{{print $2,$4,$5,$1,"0",$3,$6,$7,"0",$8,$9,$10}}' ${{bn}}.tmp > {output.bed12}
rm -f ${{bn}}.tmp
"""
