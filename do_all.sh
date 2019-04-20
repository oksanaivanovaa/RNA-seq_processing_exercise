# getting fastqs
mkdir fastqs
ln -s /mnt/GSE103958/*.fastq.gz fastqs/

TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/.fastq.gz//')

# strting with only one sample
TAG=MPH_untr_rep1

for TAG in $TAGS; do
  # fastqc
  
  OUTDIR="fastqc/$TAG"; mkdir -p "$OUTDIR"
  fastqc -o "$OUTDIR" "fastqs/$TAG.fastq.gz" |& tee "$OUTDIR/$TAG.fastqc.log"
done
  #======= Hisat2 pipeline ==========

TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)

for TAG in $TAGS; do
  HISAT_IDX=/mnt/reference/Gencode_mouse/release_M20/GRCm38.primary_assembly
  
  # aligning to the genome reference
  
  OUTDIR="hisat2/$TAG"; mkdir -p "$OUTDIR"
  date
  hisat2 -p 8 --new-summary -x  ${HISAT_IDX} \
    -1 "fastqs/$TAG*_1.fastq.gz" -2 "fastqs/$TAG*_2.fastq.gz" \
    2> "$OUTDIR/$TAG.hisat2.log" \
    | samtools view -b - > "$OUTDIR/$TAG.raw.bam"
  date
  
done  


cat "$OUTDIR/$TAG.hisat2.log"
ls $OUTDIR


  # post-processing the alignments
for TAG in $TAGS; do
  OUTDIR="hisat2/$TAG"; mkdir -p "$OUTDIR"
  date
  samtools sort -@ 8 -O bam "$OUTDIR/$TAG.raw.bam" > "$OUTDIR/$TAG.bam" && \
    samtools index "$OUTDIR/$TAG.bam" && \
    rm -v "$OUTDIR/$TAG.raw.bam"
  date
  
  # calculating coverage for vizualization
  bamCoverage -b "$OUTDIR/$TAG.bam" -o "$OUTDIR/$TAG.cov.bw" |& tee "$OUTDIR/$TAG.bamcov.log"
done



  # QC

#TAGS="SRX3195614 SRX3195615 SRX3195616 SRX3195617"

for TAG in $TAGS; do
  OUTDIR="hisat2/$TAG"
  REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
  infer_experiment.py -i "$OUTDIR/$TAG.bam" \
    -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.infer_experiment.txt"
  
  REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_Gencode_VM18.bed
  read_distribution.py -i "$OUTDIR/$TAG.bam" \
    -r $REFGENE_MODEL | tee "$OUTDIR/$TAG.read_distribution.txt"
done
  
  # REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10_rRNA.bed
  # split_bam.py -i "$OUTDIR/$TAG.bam" -r $REFGENE_MODEL -o "$OUTDIR/$TAG.rrna" | tee "$OUTDIR/$TAG.split_rrna.txt"
  # rm "$OUTDIR/$TAG.rrna.ex.bam" "$OUTDIR/$TAG.rrna.in.bam" "$OUTDIR/$TAG.rrna.junk.bam"
  
  # REFGENE_MODEL=/mnt/reference/Gencode_mouse/release_M20/mm10.HouseKeepingGenes.bed
  # geneBody_coverage.py \
  #   -i $OUTDIR/$TAG.bam \
  #   -o $OUTDIR/$TAG \
  #   -r $REFGENE_MODEL 
    



# Counting reads
TAGS=$(ls fastqs/SRX31956*.fastq.gz | xargs -n 1 basename | sed 's/_[1,2].fastq.gz//' | uniq)

## 3'libraries!!!!!
for TAG in $TAGS; do  
  GTF=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf
  
  OUTDIR="featureCounts/$TAG"; mkdir -p "$OUTDIR"
  date
  featureCounts -a "$GTF" -s 0 -p -o "$OUTDIR/$TAG.fc.txt" \
    "hisat2/$TAG/$TAG.bam" |& tee "$OUTDIR/$TAG.fc.log"
  date
done   

  head "$OUTDIR/$TAG.fc.txt"
  wc -l "$OUTDIR/$TAG.fc.txt"

#### STOP

  #========== Kallisto ======================
  
  mkdir kallisto
  
  KALLISTO_IDX=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.transcripts.kalliso.idx
  
  OUTDIR="kallisto/$TAG"; mkdir -p "$OUTDIR"
  date
  # --single -l and -s option should be set for each dataset separately, 200+-50 is most common for single end
  kallisto quant -i $KALLISTO_IDX -t 8 \
    --single -l 200 -s 50 \
    --plaintext \
    -o $OUTDIR \
    fastqs/$TAG.fastq.gz |& tee $OUTDIR/$TAG.kallisto.log
  date
  
done

#========== multiqc for everything =============

multiqc -x .Rproj.user -f .

#========== mmquant ============
OUTDIR="mmquant"; mkdir -p "$OUTDIR"
GTF=/mnt/reference/Gencode_mouse/release_M20/gencode.vM20.annotation.gtf
date
mmquant -a "$GTF" -s U -o "$OUTDIR/mmq.txt" \
  -r hisat2/*/*.bam |& tee "$OUTDIR/mmq.log"
date