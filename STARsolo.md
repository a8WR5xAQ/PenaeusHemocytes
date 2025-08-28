# STARsolo mapping
## Downloads geneome sequence
For *P. japonicus*<br>
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/017/312/705/GCF_017312705.1_Mj_TUMSAT_v1.0/GCF_017312705.1_Mj_TUMSAT_v1.0_genomic.fna.gz<br><br>
For *P. vannamei*<br>
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/789/085/GCF_003789085.1_ASM378908v1/GCF_003789085.1_ASM378908v1_genomic.fna.gz<br><br>
For *P. monodon*<br>
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/015/228/065/GCF_015228065.2_NSTDA_Pmon_1/GCF_015228065.2_NSTDA_Pmon_1_genomic.fna.gz<br><br>
For *D. melanogaster*<br>
https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz<br><br>

## Creation of STAR index
```
STAR \
--runThreadN X \
--runMode genomeGenerate \
--genomeDir Dir_to_genome \
--genomeFastaFiles file_of_genome \
--sjdbGTFfeatureExon exon \
--sjdbGTFfile GTFfile.gtf
```

## Mapping against genome sequence and count UMIs/Genes by STARsolo
```
STAR \
--runThreadN X \
--runMode alignReads \
--outSAMtype BAM SortedByCoordinate \
--sysShell /bin/bash \
--genomeDir Dir_to_index \
--readFilesIn Read2.fastq.gz Read1.fastq.gz \
--readFilesCommand 'gzip -c -d' \
--soloCBwhitelist None \
--soloType CB_UMI_Simple \
--soloCBmatchWLtype 1MM_multi \
--soloCBstart  1 \
--soloCBlen 12 \
--soloUMIstart 13 \
--soloUMIlen 8 \
--soloBarcodeReadLength 0 \
--sjdbGTFfile Dir_to_gtf \
--soloCellFilter CellRanger2.2 1500 0.99 10 \
--outFilterMultimapNmax 1 \
--outSAMattributes NH HI AS nM CB UB CR CY UR UY \
--outFileNamePrefix Dir_to_out \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--bamRemoveDuplicatesType UniqueIdentical \
--outTmpDir Dir_to_tmp \
--soloFeatures Gene
```