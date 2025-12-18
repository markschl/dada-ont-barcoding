#!/bin/bash
#
#SBATCH -A node              # GPU account
#SBATCH -p gpu               # GPU partition
#SBATCH --mem 20G            # initially a lot of memory is used, later apparently not anymore
#SBATCH -t 2-00:00           # max. 2 days
#SBATCH -n 2                 # max. 2 tasks
#SBATCH --gres=gpu:1         # request 1 GPU
#SBATCH --mail-type=END,FAIL # notify you at given events
#SBATCH --mail-user=<email>  # Your email

################################################################################
# Things to do:
# - make sure that the above settings are correct for your system
# - make sure that Dorado is installed
# - adjust the Dorado path below
# - add your email address for notifications (above: --mail-user=<email>)
################################################################################

set -euo pipefail

nvidia-smi >&2

model="$1"
indir="$2"

dorado=~/opt/dorado/bin/dorado  # <- adjust to install location

"$dorado" basecaller $model --verbose --no-trim --emit-fastq "$indir" | 
  /bin/gzip -c > output.fastq.gz
