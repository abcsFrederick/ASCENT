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


rule rmats:
    input:
        unpack(get_samples_files)
    output:
        summary=join(RESULTSDIR,"{contrast}","rmats","summary.txt")
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
