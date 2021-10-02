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