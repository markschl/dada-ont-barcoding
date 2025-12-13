
# Basecalling

If the computer has a Nvidia GPU with Cuda support (or a connected [eGPU](https://egpu.io)), then basecalling may be done locally.
Alternatively, basecalling may be done on a different computer, e.g. a HPC server with a GPU node.

## Model

Different [machine learning models](https://software-docs.nanoporetech.com/dorado/latest/models/models) are available
for the basecalling, offering different trade-offs between speed and accuracy.
Our experience suggests that the `sup` model gives substantially better results than the `hac` model (by far) in terms of Q scores
and recognized primers/barcodes. Basecalling can sometimes take a full day for a MinION run, depending on the GPU setup.

## On your computer

[Dorado](https://software-docs.nanoporetech.com/dorado/latest) is needed for this (`dorado` command should be available in console). It's installation size is quite large (a few Gigabytes).

If not using Dorado for other things, you may also use the `install-dorado.sh` installation script, which automatically downloads the software into the *bin* subfolder.
Make sure to supply the latest version, which can be [checked on the website](https://software-docs.nanoporetech.com/dorado/latest) (example with `1.3.0`):

```sh
cd path/to/dada-ont-barcoding
scripts/install-dorado.sh 1.3.0
```

The basecalling is done as follows:

```sh
cd path/to/dada-ont-barcoding

# create the input directory (if not present)
run_name=my_run
mkdir -p analysis/$run_name

# Directory with the Nanopore data (may be on an external hard drive)
ont_data=/media/user/run/data

# The basecalling model
model=sup

# run the base-caller (may take a long time)
# 'bin/dorado' if Dorado was installed with 'scripts/install-dorado.sh'
dorado basecaller $model -r --no-trim --emit-fastq $ont_data | gzip -c > analysis/$run_name/reads.fastq.gz
```

An easy way to find out the actual paths of the nanopore data is to drag & drop the folders on the console, then copy the path there.


### On a HPC server

Basecalling on a remote server is slightly more complicated. This script shows how to do this on a SLURM server with a GPU node. The `remote_dir` must contain a job submission script (based on this [template](https://github.com/markschl/dada-ont-barcoding/templates/base_calling/slurm_submit_basecaller.sh)); it can be copied to the remote directory e.g. using an FTP client; adjust the `dorado` path and your email address inside the script.

```sh
cd path/to/ont-barcoding

# run directory
run_name=my_run
run_dir=analysis/$run_name

# SSH login to HPC server
ssh=user@server-address

# Remote directory containing the submission script (see templates/slurm_submit_basecaller.sh)
# The base-called data is placed in $remote_dir/$run_name
remote_dir=projects/ONT

# Nanopore data export directory
ont_data=/media/user/run/data

# The basecalling model
model=sup

# copy the data
cd $run_dir
rsync -av --info=progress2 "$ont_data"/*/*/pod5_skip $ssh:$remote_dir/$run_name
# submit the run
ssh $ssh -t "cd ~/$remote_dir && sbatch slurm_submit_basecaller.sh $model $run_name/pod5_skip"
```

There should be a message: "submitted batch job ...". Wait for the completion email and then execute:

```sh
rsync -av $ssh/$remote_dir/$run_name/output.fastq.gz $analysis/$run_name/reads.fastq.gz
```
