#!/bin/bash
#SBATCH --job-name=MEIsensor_demo
#SBATCH --output=demo/result/demo_run.%j.log
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --time=00:30:00

source ~/anaconda3/etc/profile.d/conda.sh
conda activate MEIsensor

# annotation (--annotate): set SVAN_PATH and SVAN_RESOURCES to enable; see README "Annotation"
REF=/path/to/GRCh38.fasta
MODEL=src/trained_model/seq_classify_model.ph
BAM=demo/data/demo.bam

SVAN_PATH=${SVAN_PATH:-}
SVAN_RESOURCES=${SVAN_RESOURCES:-}

OUTDIR=demo/result
mkdir -p $OUTDIR

ANNOT_ARGS=""
if [ -n "$SVAN_PATH" ] && [ -n "$SVAN_RESOURCES" ]; then
    ANNOT_ARGS="--annotate true --svan-path $SVAN_PATH --svan-resources $SVAN_RESOURCES"
fi

srun python src/main.py \
    -i $BAM \
    -v $OUTDIR/demo.MEI.vcf \
    --reference $REF \
    -m $MODEL \
    -t ${SLURM_CPUS_PER_TASK:-8} \
    $ANNOT_ARGS
