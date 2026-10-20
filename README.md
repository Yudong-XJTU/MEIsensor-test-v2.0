# MEIsensor (v1.0)

MEIsensor is a tool for detecting and classifying mobile element insertions (MEIs) from long-read sequencing data. It calls insertions from aligned reads and uses a trained convolutional neural network to classify each insertion as Alu, LINE1 or SVA.

## Overview

MEIsensor extracts split-read and inline insertion signals from a coordinate-sorted BAM, clusters them into candidate insertions, and assembles a consensus insertion sequence from the supporting reads. Each candidate is then classified by the neural network; insertions that are not classified as MEI are dropped, so the output VCF only contains Alu/LINE1/SVA calls.

![Overview of MEIsensor](figure/figure1.png)

## Installation

Tested on Linux with Python 3.10, PyTorch 2.7 (CUDA 12.8) and an NVIDIA GPU.

```bash
git clone https://github.com/Yudong-XJTU/MEIsensor-v1.0.git
cd MEIsensor-v1.0
```

With CUDA 12.8 or higher:

```bash
conda env create -f environment.yml -n MEIsensor
conda activate MEIsensor
```

For other CUDA versions, install manually and pick the PyTorch build that matches yours:

```bash
conda create -n MEIsensor python=3.10
conda activate MEIsensor
pip install torch==2.7.0 torchvision==0.22.0 torchaudio==2.7.0 --index-url https://download.pytorch.org/whl/cu128
pip install -r requirements.txt
```

Check CUDA:

```bash
python -c "import torch; print(torch.__version__, torch.version.cuda, torch.cuda.is_available())"
```

## Usage

```bash
python src/main.py -i $BAM_PATH -v $OUTPUT_PATH --reference $REF_PATH -m $MODEL_PATH -t $THREAD
```

Example:

```bash
python src/main.py \
    -i sample.hifi.minimap2.sorted.bam \
    -v sample.MEI.vcf \
    --reference GRCh38.fasta \
    -m src/trained_model/seq_classify_model.ph \
    -t 12
```

Parameters:

```text
-i BAM_PATH                       coordinate-sorted, indexed .bam of aligned long reads
-v OUTPUT_PATH                    output VCF (.vcf or bgzipped .vcf.gz)
--reference REF_PATH              reference genome the reads were aligned to
-m MODEL_PATH                     trained model, shipped at src/trained_model/seq_classify_model.ph
-t THREAD                         number of threads (default 4)
```

The reference is needed for classification: the base at each insertion site is taken from it and prepended to the ALT sequence before the model is applied. Without a reference no classified calls are written. Use the same reference the BAM was aligned against.

## Output

The output is a standard VCF:

- ID: `MEIsensor.ALU.INS`, `MEIsensor.LINE1.INS` or `MEIsensor.SVA.INS`
- REF/ALT: reference base at the insertion point, and the reference base plus the assembled insertion sequence
## Demo

`demo/data/` contains a small BAM extracted from HG00514 (chr1:28-29Mb). Run it from the repository root:

```bash
python src/main.py -i demo/data/HG00514.chr1_28-29Mb.demo.bam -v demo.MEI.vcf \
    --reference GRCh38.fasta -m src/trained_model/seq_classify_model.ph -t 8
```

This takes a few seconds and calls 4 Alu insertions. `demo/result/` holds the expected output.

On a SLURM cluster, set `REF` in `demo/run_demo.sbatch` to your reference genome, then:

```bash
sbatch demo/run_demo.sbatch
```
