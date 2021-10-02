rule create_star_index:
    input:
        reffa=REFFA,
        gtf=GTF
    output:
        sa=join(STARINDEXDIR,"SA")
    params:
        rl=MAXRL,
        stardir=STARINDEXDIR
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
    --outTmpDir /lscratch/${{SLURM_JOB_ID}}/tmp_{params.rl}
"""

rule star:
    input:
        R1=join(WORKDIR,"fastqs","{sample}.R1.fastq.gz"),
        R2=join(WORKDIR,"fastqs","{sample}.R2.fastq.gz"),
        sa=rules.create_star_index.output.sa
    output:
        junction=join(BAMDIR,"{sample}.Chimeric.out.junction"),
        bam=join(BAMDIR,"{sample}.Aligned.sortedByCoord.out.bam"),
        tbam=join(BAMDIR,"{sample}.Aligned.toTranscriptome.out.bam"),
        genecounts=join(BAMDIR,"{sample}.ReadsPerGene.out.tab"),
    threads: getthreads("star")
    params:
        sample="{sample}",
        gtf=GTF,
        rl=MAXRL,
        starindexdir=STARINDEXDIR,
        bamdir=join(WORKDIR,"bam")
    envmodules: 
        TOOLS['rmats']['version']
    shell:"""
set -euf -o pipefail
overhang=$(echo {params.rl}|awk '{{print $1-1}}')
cd {params.bamdir}
STAR \
    --chimSegmentMin 2 \
    --outFilterMismatchNmax 3 \
    --alignEndsType EndToEnd \
    --runThreadN {threads} \
    --outSAMstrandField intronMotif \
    --outSAMtype BAM SortedByCoordinate \
    --alignSJDBoverhangMin 6 \
    --alignIntronMax 299999 \
    --genomeDir {params.starindexdir} \
    --sjdbGTFfile {params.gtf} \
    --outFileNamePrefix "{params.sample}." \
    --readFilesIn {input.R1} {input.R2} \
    --readFilesCommand zcat \
    --outTmpDir=/lscratch/${{SLURM_JOB_ID}}/{params.sample} \
    --sjdbOverhang $overhang \
    --quantMode GeneCounts TranscriptomeSAM
"""
