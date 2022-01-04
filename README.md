# CCBR_ASE
ASE (alternate splicing events) identification and quantification workflow for multi-group multi-contrasts scenarios.
Runs:

  * rMATS
  * IsoformSwitchAnalyzeR

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

## Note:
IsoformSwitchAnalyzeR fails if replicates are absent.
