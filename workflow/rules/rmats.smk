def get_samples_files(wildcards):
    files=dict()
    contrast=wildcards.contrast
    x=CONTRASTSDF[CONTRASTSDF['name']==contrast]
    x=x.iloc[0]
    g1=x.group1
    files['s1']=GROUP2SAMPLESTXT[g1]
    g2=x.group2
    files['s2']=GROUP2SAMPLESTXT[g2]
    return files

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

rule rmats:
    input:
        unpack(get_samples_files),
	sa=rules.create_star_index.output.sa
    output:
        summary=join(RESULTSDIR,"{contrast}","summary.txt")
    params:
        rl=MAXRL,
        gtf=GTF,
        stardir=STARINDEXDIR,
    envmodules:
        TOOLS['rmats']['version']
    threads: getthreads("rmats")
    shell:"""
set -euf -o pipefail
outdir=$(dirname {output.summary})
cp "{input.s1}" $outdir/
cp "{input.s2}" $outdir/
python ${{RMATS_SRC}}/rmats.py \
    --s1 "{input.s1}" \
    --s2 "{input.s2}" \
    --od "$outdir" \
    --gtf "{params.gtf}" \
    --readLength {params.rl} \
    --bi "{params.stardir}" \
    --nthread {threads} \
    -t "paired" \
    --novelSS \
    --libType "fr-secondstrand" \
    --tmp /lscratch/${{SLURM_JOB_ID}}/
"""
