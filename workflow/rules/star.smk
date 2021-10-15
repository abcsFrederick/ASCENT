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

rule get_rsem_counts:
    input:
        bed12=rules.create_bed12.output.bed12,
        bam=rules.star.output.bam,
        tbam=rules.star.output.tbam        
    output:
        strandinfo=join(WORKDIR,"strandinfo","{sample}.strandinfo"),
        gcounts=join(WORKDIR,"rsem","genecounts","{sample}","{sample}.RSEM.genes.results"),
        tcounts=join(WORKDIR,"rsem","isoformcounts","{sample}","{sample}.RSEM.isoforms.results"),
    envmodules: TOOLS['rseqc']['version'], TOOLS['rsem']['version'], TOOLS["samtools"]["version"]
    threads: getthreads("get_rsem_counts")
    params:
        sample="{sample}",
        rsemdir=join(WORKDIR,"rsem")
    shell:"""
set -euf -o pipefail
cd {params.rsemdir}
# inter strandedness
samtools index -@{threads} {input.bam}
infer_experiment.py -r {input.bed12} -i {input.bam} -s 1000000 > {output.strandinfo}
# Get strandedness to calculate Forward Probability
fp=$(tail -n1 {output.strandinfo} | awk '{{if($NF > 0.75) print "0.0"; else if ($NF<0.25) print "1.0"; else print "0.5";}}')
echo "Forward Probability Passed to RSEM: $fp"
rsemindex=$(echo {input.bed12}|sed "s@.bed12@@g")
rsem-calculate-expression --no-bam-output --calc-ci --seed 12345  \
        --bam --paired-end -p {threads}  {input.tbam} $rsemindex {params.sample} --time \
        --temporary-folder /lscratch/$SLURM_JOBID --keep-intermediate-files --forward-prob=${{fp}} --estimate-rspd
mv {params.sample}.genes.results {output.gcounts}
mv {params.sample}.isoforms.results {output.tcounts}
"""    

