rule create_star_index:
    input:
        reffa=REFFA,
        gtf=GTF
    output:
        sa=join(STARINDEXDIR,"SA")
    params:
        rl=MAXRL,
        stardir=STARINDEXDIR,
    envmodules: 
        TOOLS['rmats']['version']
    threads: getthreads("create_star_index")
    shell:"""
set -euf -o pipefail
mkdir -p {params.stardir}
STAR \
    --runThreadN {threads} \
    --runMode genomeGenerate \
    --genomeDir {params.stardir} \
    --genomeFastaFiles {input.reffa} \
    --sjdbGTFfile {input.gtf} \
    --sjdbOverhang $(echo {params.rl}|awk '{{print $1-1}}') \
    --outTmpDir /lscratch/${{SLURM_JOB_ID}}/tmp_{params.rl}
"""
