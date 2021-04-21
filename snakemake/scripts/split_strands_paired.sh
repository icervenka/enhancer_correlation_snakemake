#!/usr/bin/env sh
# Adapted from Istvan Albert - https://www.biostars.org/p/92935/
input=$(realpath $1)
outdir=$(realpath $2)
outputf=$(realpath $3)
outputr=$(realpath $4)
ebed=$(realpath $5)
threads=$6

# Forward strand.
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
samtools view -@ $threads -L $ebed -b -f 128 -F 16 $input > $outdir/fwd1.temp.bam
samtools index $outdir/fwd1.temp.bam

samtools view -@ $threads -L $ebed -b -f 80 $input > $outdir/fwd2.temp.bam
samtools index $outdir/fwd2.temp.bam

# Combine alignments that originate on the forward strand.
samtools merge -@ $threads -f $outputf $outdir/fwd1.temp.bam $outdir/fwd1.temp.bam
samtools index $outputf

# Reverse strand
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
samtools view -@ $threads -L $ebed -b -f 144 $input > $outdir/rev1.temp.bam
samtools index $outdir/rev1.temp.bam

samtools view -@ $threads -L $ebed -b -f 64 -F 16 $input > $outdir/rev2.temp.bam
samtools index $outdir/rev2.temp.bam

# Combine alignments that originate on the reverse strand.
samtools merge -@ $threads -f $outputr $outdir/rev1.temp.bam $outdir/rev2.temp.bam
samtools index $outputr

rm -f $outdir/fwd1.temp.bam
rm -f $outdir/fwd1.temp.bam.bai
rm -f $outdir/fwd2.temp.bam
rm -f $outdir/fwd2.temp.bam.bai

rm -f $outdir/rev1.temp.bam
rm -f $outdir/rev1.temp.bam.bai
rm -f $outdir/rev2.temp.bam
rm -f $outdir/rev2.temp.bam.bai
