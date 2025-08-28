# Kallisto mapping
## Creation of Kalisto index
```
kallisto index -i path/to/index path/to/genome
```

## Count Genes by Kallisto
```
fastq='path/to/fastq'
out='path/to/outdir'
for id in `cut -f 1 -d ',' path/to/SRA_list.csv`
do
    mkdir ${fastq}/${id}
    cd ${fastq}/${id}
    wget https://sra-pub-run-odp.s3.amazonaws.com/sra/${id}/${id}
    fasterq-dump ${id}

    mkdir ${out}/${id}
    kallisto quant -i path/to/index -o ${out}/${id} ${fastq}/${id}/${id}_1.fastq ${fastq}/${id}/${id}_2.fastq
    rm ${fastq}/${id}/*
done
```