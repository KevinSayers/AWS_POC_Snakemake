1. Create a fastq folder and add the compressed RNA-seq fastq files to the directory. The files must contain the extension _R1.fastq.gz & _R2.fastq.gz
2. The wrapper run_SN_MVP.sh is used to launch the pipeline and to do a dry run. 
```
% ./run_SN_MVP.sh
Pipeline Dir: /vf/users/khanlab/projects/AWS_MVP_HPC/snakemake
Git Commit/Tag: 20675a530d7de78f548fbba0c71b1ad4d8950369	

run_SN_MVP.sh
--> run Snakemake_MVP pipeline

USAGE:
  bash ./run_SN_MVP.sh -m/--runmode=<RUNMODE> -w/--workdir=<WORKDIR>
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

**init**  Initialize the output folder

Specify an output folder and initialize it with the following command. 

```
bash ./run_SN_MVP.sh -m=init -w=<path to Output Directory>
```
**dryrun** To verify the run

```cd``` to the output directory and edit the config.yaml. Add the run specific user defined data like the input files location, path to references and genome annotation, fasta files. Then we run the following command to ensure that we are ready to submit the job to the cluster.

``` 
bash run_SN_MVP.sh -m=dryrun -w=<path to Output Directory>
```
**run** launching the pipeline

If the dryrun looks good, we can go ahead and launch the pipeline using the following command.

```
bash run_SN_MVP.sh -m=run -w=<path to Output Directory>
```
