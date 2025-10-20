# Characterization of the dynamic microbiome evolution across thrips species
A pipeline for analyzing the microbiome of thrips

## characterization of thrips microbiome

### short reads

Kraken2

SingleM


### long reads
To assemble the bacterial genomes from long reads, we first removed the host reads:

### a. host read removal

    
```
conda activate samtools_v1.21
# for HiFi reads
minimap2 -a genome.fa -o genome.sam -t 20 -x map-hifi hifi.fq.gz

# for ONT reads
minimap2 -a genome.fa -o genome.sam -t 20 -x map-ont ont.fq.gz

# for CLR reads
minimap2 -a genome.fa -o genome.sam -t 20 -x map-pb clr.fq.gz

samtools view -h -b -S genome.sam -o genome.bam --threads 20
samtools sort genome.bam -o genome.sorted.bam --threads 20
samtools index genome.sorted.bam
```

myloasm

metaflye



## genome assembly of bacterial isolates
