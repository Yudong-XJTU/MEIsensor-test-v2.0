#!/bin/bash
#SBATCH --job-name=MEIsensor_demo
#SBATCH --output=demo/result/demo_run.%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00

source ~/anaconda3/etc/profile.d/conda.sh
conda activate MEIsensor

REF=/path/to/GRCh38.fasta
MODEL=src/trained_model/seq_classify_model.ph
BAM=demo/data/HG00514.chr1_28-29Mb.demo.bam

OUTDIR=demo/result
mkdir -p $OUTDIR

srun python src/main.py \
    -i $BAM \
    -v $OUTDIR/demo.MEI.vcf \
    --reference $REF \
    -m $MODEL \
    -t ${SLURM_CPUS_PER_TASK:-8}
