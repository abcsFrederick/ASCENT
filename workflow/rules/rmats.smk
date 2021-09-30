def get_fastq_R1s():
    R1s=dict()
    for sample in SAMPLESDF.index:
        R1s[sample]=SAMPLESDF["R1"][sample]
    return(R1s)


rule get_max_rl:
    input: 
        unpack(get_fastq_R1s)
    output: 
        maxrl=join(RESULTSDIR,"maxrl.txt")
    shell:"""
set -euf -o pipefail
max=0
for f in {input};do
    l2=$(zcat $f|head -n2|tail -n1|wc -c)
    l2=$(echo $l2|awk '{{print $1-1}}')
    if [ "$l2" -gt "$max" ];then
        max=$l2
    fi
done
echo "$max" > {output}
"""

rule create_star_index:
    input:
        reffa=REFFA,
        gtf=GTF
    output:
        sa=join(STARINDEXDIR,"SA")
    params:
        rl=MAXRL,
        stardir=STARINDEXDIR,
    envmodules: TOOLS['rmats']['version']
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
    --sjdbOverhang $(echo "{params.rl}"|awk '{{print $1-1}}') \
    --outTmpDir /lscratch/${{SLURM_JOB_ID}}/tmp_{params.rl}
"""