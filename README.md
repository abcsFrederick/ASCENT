# ASCENT

**A**lternate **S**pli**C**ing **E**ve**N**t **T**ools

[![release](https://img.shields.io/github/v/release/CCBR/ASCENT?color=blue&label=latest%20release)](https://github.com/CCBR/ASCENT/releases/latest)

ASE (alternate splicing events) are identified and quantified using the ASCENT (**A**lternate **S**pli**C**ing **E**ve**N**t **T**ools) pipeline. This workflow can be used for multi-group multi-contrasts scenarios.

ASCENT currently runs:

  * [rMATS<sup>1</sup>](http://rnaseq-mats.sourceforge.net/)
  * [IsoformSwitchAnalyzeR<sup>2</sup>](https://bioconductor.org/packages/release/bioc/html/IsoformSwitchAnalyzeR.html)

> **NOTE**: 
> More ASE tools will be added to ASCENT in the future as and when needed. Please contact [Vishal Koparde](mailto:vishal.koparde@nih.gov) to make requests to add specific ASE tool(s).

## Usage:

```
 % bash /data/Ziegelbauer_lab/Pipelines/CCBR_ASE/v0.1.1/run
Pipeline Dir: /gpfs/gsfs12/users/Ziegelbauer_lab/Pipelines/CCBR_ASE/v0.1.1
Git Commit/Tag: a41b33007b6063b1d95744808f527a126fb1f0dc	v0.1.1
/gpfs/gsfs12/users/Ziegelbauer_lab/Pipelines/CCBR_ASE/v0.1.1/run
--> run CCBR Alternate Splicing Events Pipeline

USAGE:
  bash /data/Ziegelbauer_lab/Pipelines/CCBR_ASE/v0.1.1/run -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
Required Arguments:
1.  RUNMODE: [Type: String] Valid options:
    *) init : initialize workdir
    *) run : run with slurm
    *) reset : DELETE workdir dir and re-init it
    *) dryrun : dry run snakemake to generate DAG
    *) unlock : unlock workdir if locked by snakemake
    *) runlocal : run without submitting to sbatch
2.  WORKDIR: [Type: String]: Absolute or relative path to the output folder with write permissions.
```

> ## **DISCLAIMER**:
> IsoformSwitchAnalyzeR fails if replicates are absent.

## References:

- <sup>1</sup> [http://dx.doi.org/10.1073/pnas.1419161111](http://dx.doi.org/10.1073/pnas.1419161111)
- <sup>2</sup> [https://doi.org/10.1158/1541-7786.MCR-16-0459](https://doi.org/10.1158/1541-7786.MCR-16-0459)


