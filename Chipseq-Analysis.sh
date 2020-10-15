#usr/bin/bash
##############
# the code was used to analyze Chip-Seq and Cut&run data
##############

#### Quality Control
ls *.gz | while read id; do fastqc $id; done

#### Trimmomatic
ls *_1.fq.gz | while read file; name=$(basename $file | cut -d _ -f 1);java -jar trimmomatic-0.39.jar PE -threads 8 $file ${name}_2.fastq Trimmomatic/${name}_1.p.fastq Trimmomatic/${name}_1.u.fq Trimmomatic/${name}_2.p.fastq Trimmomatic/${name}_2.u.fq ILLUMINACLIP:adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:4:20 MINLEN:40 HEADCROP:10; done

#### bowtie2
ls *_1.p.fastq | while read id; do name=$(basename $id | cut -d _ -f 1); bowtie2 -x mm10 -1 $id -2 ${name}_2.p.fastq -t -q -N 1 -L 25 -S ${name}.sam; done

#### samtools
ls *.sam | while read id; do name=$(basename $id | cut -d . -f 1); perl get_bt2_sam_uniq.pl $id ${name}.uniq.sam 30; done
ls *.uniq.sam | while read id; do name=$(basename $id | cut -d . -f 1); samtools view -h -F 4 -bS -o ${name}.bam -@ 3 $id; done
ls *.bam | while read id; do name=$(basename $id | cut -d . -f 1); samtools sort -@ 4 -o ${name}.sorted.bam $id; done
ls *.sorted.rmdup.bam | while read id; do samtools index $id ; done
ls *.sorted.bam | while read id; do name=$(basename $id | cut -d . -f 1); java -jar ~/Software/picard/build/libs/picard.jar MarkDuplicates I=$id O=${name}.sorted.rmdup.bam REMOVE_DUPLICATES= true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT  M=${name}.sort.addhead.rmdup.metric; done   


#### SNPSplit
cat file.txt | while read id; do echo $id; bowtie2 -t -q -p 4 -N 1 -L 25 -X 2000 --no-mixed --no-discordant -x N-mask-genome/PWK_PhJ_N-masked/PWK_PhJ_N_Genome -1 $id/${id}_1.clean.fq -2 $id/${id}_2.clean.fq -S ${id}.sam; done

ls *.uniq.sam | while read id; do name=$(basename $id | cut -d . -f 1); samtools view -h -F 4 -bS -o ${name}.bam -@ 3 $id; done
ls *.bam | while read id; do name=$(basename $id | cut -d . -f 1); samtools sort -@ 4 -o ${name}.sorted.bam $id; done
ls *.sorted.rmdup.bam | while read id; do samtools index $id ; done
ls *.sorted.bam | while read id; do name=$(basename $id | cut -d . -f 1); java -jar ~/Software/picard/build/libs/picard.jar MarkDuplicates I=$id O=${name}.sorted.rmdup.bam REMOVE_DUPLICATES= true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000 VALIDATION_STRINGENCY=LENIENT  M=${name}.sort.addhead.rmdup.metric; done

ls PM/*.sorted.bam | while read id; do SNPsplit --paired --snp_file snp_file $id --confilcting; done

#### Bedtools
ls *.sorted.rmdup.genome2.bam | while read id; do name=$(basename $id | cut -d . -f 1); bamToBed -i $id > Bedtools/${name}.genome2.bed; done

ls *.sorted.rmdup.genome1.bam | while read id; do name=$(basename $id | cut -d . -f 1); bamToBed -i $id > Bedtools/${name}.genome1.bed; done

ls *.sorted.rmdup.bam | while read id; do name=$(basename $id | cut -d . -f 1); bamToBed -i $id > Bedtools/${name}.bed; done

#### Deeptools
plotFingerprint -b LXY-0109-8_FKDL192507604-1a_HYG3GCCXY_L6.sorted.rmdup.bam LXY-0109-9_FKDL191548618-1a_HYG3GCCXY_L6.sorted.rmdup.bam --labels H3K27m2 input --minMappingQuality 30 --skipZeros --numberOfSamples 50000 -T "WT Fingerprints" --plotFile WT-08-09-fingerprint.pdf --outRawCounts WT-08-09-fingerprint.tab

plotFingerprint -b LXY-0109-17_FKDL191549261-1a_HYG3GCCXY_L6.sorted.rmdup.bam LXY-0109-18_FKDL191549262-1a_HYG3GCCXY_L6.sorted.rmdup.bam --labels H3K27m2 input --minMappingQuality 30 --skipZeros --numberOfSamples 50000 -T "WT Fingerprints" --plotFile WTZygote-17-18-fingerprint.pdf --outRawCounts WTZygote-17-18-fingerprint.tab

plotCoverage -b *.sorted.rmdup.bam --plotFile WT-coverage.pdf --numberOfSamples 500000  --outRawCounts coverage.tab -T “Coverage by Deeptools” –skipZeros 

ls *sorted.rmdup.bam | while read id; do name=$(basename $id | cut -d . -f 1); bamCoverage -b $id --normalizeUsing RPKM --binSize 30 --smoothLength 300 -p 4 --extendReads 200 -o Deeptools/${name}.bw; done

ls *sorted.rmdup.bam | while read id; do name=$(basename $id | cut -d . -f 1); bamCoverage -b $id -o ${name}.bw; done

ls *sorted.rmdup.bam | while read id; do name=$(basename $id | cut -d . -f 1); bamCoverage -e 170 -bs 100 -of bedgraph -b  -o ${name}.bdg; done


#### SICER
/home/mengqingren/miniconda3/envs/sicer/bin/SICER-df.sh ../Bedtools/LXY-1130-18_TKD181200076_HVM7CCCXY_L3.bed ../Bedtools/LXY-1130-19_TKD181200077_HVM7CCCXY_L3.bed ../Bedtools/LXY-0109-8_FKDL192507604-1a_HYG3GCCXY_L6.bed ../Bedtools/LXY-0109-9_FKDL191548618-1a_HYG3GCCXY_L6.bed 200 600 1e-5 1e-5

/home/mengqingren/miniconda3/envs/sicer/bin/SICER-df.sh ../Bedtools/LXY-1130-18_TKD181200076_HVM7CCCXY_L3.bed ../Bedtools/LXY-1130-19_TKD181200077_HVM7CCCXY_L3.bed ../Bedtools/LXY-0109-8_FKDL192507604-1a_HYG3GCCXY_L6.bed ../Bedtools/LXY-0109-9_FKDL191548618-1a_HYG3GCCXY_L6.bed 200 600 1e-5 1e-5

SICER.sh Bedtools LXY-1130-18_TKD181200076_HVM7CCCXY_L3.bed  LXY-1130-19_TKD181200077_HVM7CCCXY_L3.bed M2EZH2-200-1e-5/ mm10 1 200 150 0.74 600 1e-5

#### the peak data were used to be annotaed in R