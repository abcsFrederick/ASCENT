def get_samples_files(wildcards):
    files=dict()
    contrast=wildcards.contrast
    x=CONTRASTSDF[CONTRASTSDF['name']==contrast]
    x=x.iloc[0]
    g1=x.group1
    g2=x.group2
    for group in [g1,g2]:
        for sample in GROUP2SAMPLES[group]:
            bamfile=join(BAMDIR,sample+".Aligned.sortedByCoord.out.bam")
            files[sample]=bamfile
    files['b1']=GROUP2RMATSTXT[g1]
    files['b2']=GROUP2RMATSTXT[g2]
    files['sa']=join(STARINDEXDIR,"SA")
    return files

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



rule rmats:
    input:
        unpack(get_samples_files)
    output:
        summary=join(RESULTSDIR,"{contrast}","summary.txt")
    threads: getthreads("rmats")
    params:
        rl=MAXRL,
        gtf=GTF,
        starindexdir=STARINDEXDIR
    envmodules: 
        TOOLS['rmats']['version']
    shell:"""
set -euf -o pipefail
outdir=$(dirname {output.summary})
cp "{input.b1}" $outdir/
cp "{input.b2}" $outdir/
python ${{RMATS_SRC}}/rmats.py \
    --b1 "{input.b1}" \
    --b2 "{input.b2}" \
    --od "$outdir" \
    --gtf "{params.gtf}" \
    --readLength {params.rl} \
    --bi "{params.starindexdir}" \
    --nthread {threads} \
    -t "paired" \
    --novelSS \
    --libType "fr-secondstrand" \
    --tmp /lscratch/${{SLURM_JOB_ID}}/
"""
