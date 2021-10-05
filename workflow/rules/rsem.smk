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
genesgtf={input.GTF}
bn=$(basename $genesgtf)
gtfToGenePred -genePredExt -geneNameAsName2 $genesgtf ${{bn}}.tmp
awk '{print $2"\t"$4"\t"$5"\t"$1"\t0\t"$3"\t"$6"\t"$7"\t0\t"$8"\t"$9"\t"$10}' ${{bn}}.tmp > {output.bed12}
rm -f ${{bn}}.tmp
"""
