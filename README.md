# Characterization of the dynamic microbiome evolution across thrips species
An analytical pipeline to survey the microbiome of thrips. 

<img width="1143" height="576" alt="image" src="https://github.com/user-attachments/assets/0deecc13-68af-4ac5-8e15-4cea1b8cdfcc" />



## 1 - Characterization of thrips microbiome
To understand the thrips microbiome, we first identify bacterial species based on short reads.  
### 1.1 - Short-reads-based approach

#### a - Quality control of short reads
We first downloaded whole-genome sequencing data for thrips. Then, we removed low-quality reads using FASTP v0.23.4 https://github.com/OpenGene/fastp. 
```
# Short reads were preprocessed and quality-controlled using fastp, an ultrafast all-in-one tool. https://github.com/OpenGene/fastp
# download specified version, i.e. fastp v0.23.4

wget http://opengene.org/fastp/fastp.0.23.4
mv fastp.0.23.4 fastp
chmod a+x ./fastp
fastp -i sample1_1.fastq.gz -I sample1_2.fastq.gz -o sample1_1.clean.fastq.gz -O sample1_2.clean.fastq.gz
fastp -i sample2_1.fastq.gz -I sample2_2.fastq.gz -o sample2_1.clean.fastq.gz -O sample2_2.clean.fastq.gz
fastp -i sample3_1.fastq.gz -I sample3_2.fastq.gz -o sample3_1.clean.fastq.gz -O sample3_2.clean.fastq.gz
```

#### b - Species assignment using Kraken2
To classify the short reads, we first classified short reads based on kmer matches using Kraken2.  
```
mamba activate kraken2
time kraken2  --threads 52 --db nt --paired sample1_1.clean.fastq.gz sample1_2.clean.fastq.gz --report sample1.kreport >sample1.kraken
time kraken2  --threads 52 --db nt --paired sample2_1.clean.fastq.gz sample2_2.clean.fastq.gz --report sample2.kreport >sample2.kraken
time kraken2  --threads 52 --db nt --paired sample3_1.clean.fastq.gz sample3_2.clean.fastq.gz --report sample3.kreport >sample3.kraken

# KrakenTools is a suite of scripts to be used alongside the Kraken, KrakenUniq, Kraken 2, or Bracken programs. These scripts are designed to help Kraken users with downstream analysis of Kraken results. https://github.com/jenniferlu717/KrakenTools/
# kreport2mpa.py This program takes a Kraken report file and prints out a mpa (MetaPhlAn) -style TEXT file

python kreport2mpa.py -r sample1.kreport -o sample1.MPA.txt
python kreport2mpa.py -r sample2.kreport -o sample2.MPA.txt
python kreport2mpa.py -r sample3.kreport -o sample3.MPA.txt

#combine_mpa.py This program combines multiple outputs from kreport2mpa.py. Files to be combined must have been generated using the same kreport2mpa.py options.

python combine_mpa.py -i sample1.MPA.txt sample2.MPA.txt sample3.MPA.txt  -o 3samples.txt

# kraken2alpha.R https://github.com/YongxinLiu/EasyMicrobiome/blob/master/script/

Rscript kraken2alpha.R -i 3samples.txt

# Microbiome Helper An assortment of scripts to help process and automate various microbiome and metagenomic bioinformatic tools.https://github.com/LangilleLab/microbiome_helper/tree/master 

perl metaphlan_to_stamp.pl 3samples.norm.txt >3samples.norm.txt.spf
Rscript metaphlan_hclust_heatmap.R -i 3samples.norm.txt.spf -t Species -n 25 -o heatmap_Species
```

#### c - Species assignment using SingleM 
We then assigned reads based on read matches to highly conserved regions of Bacteria using SingleM (https://github.com/wwood/singlem)
```
conda activate singlem_v0.19.0
singlem pipe -1 read.1.fq.gz -2 read.2.fq.gz -p output --threads 20
```

### 1.2 - Long-reads-based approach
To assemble the bacterial genomes from long reads, we first removed the host reads. 

#### a - Host read removal
To remove reads from the host, we first map long reads to host genomes. Only chromosome-level scaffolds were used for mapping here to avoid mis-aligned bacterial reads. We used three different parameters in Minimap2 to align (1) PacBio HiFi reads, (2) Nanopore reads, and (3) PacBio CLR reads. 
    
```
conda activate samtools_v1.21
# 1 - for HiFi reads
minimap2 -a genome.fa -o genome.sam -t 20 -x map-hifi hifi.fq.gz

# 2 - for ONT reads
minimap2 -a genome.fa -o genome.sam -t 20 -x map-ont ont.fq.gz

# 3 - for CLR reads
minimap2 -a genome.fa -o genome.sam -t 20 -x map-pb clr.fq.gz

samtools view -h -b -S genome.sam -o genome.bam --threads 20
samtools sort genome.bam -o genome.sorted.bam --threads 20
samtools index genome.sorted.bam
```
#### b - Metagenomic assembly
To assemble bacterial genomes, we used two strategies. For PacBio HiFi reads and Nanopore reads, myloasm (https://github.com/bluenote-1577/myloasm) was used. For PacBio CLR reads, metaflye (https://github.com/mikolmogorov/Flye) was used. 

```
# 1 - PacBio HiFi reads and Nanopore reads
conda activate myloasm_v0.2.0
myloasm bac_reads.fastq.gz -o asm -t 96 --hifi

# 2 - PacBio CLR reads
conda activate flye_v2.9.6
flye --pacbio-raw bac_reads.fastq.gz --out-dir step1_flye_out --threads 96 --meta
```

#### c - Binning
To have an accurate estimation of bacterial diversity and richness, we binned metagenomic assemblies using SemiBin (https://github.com/BigDataBiology/SemiBin). 
```
# 1 - we first aligned reads to the metagenomic scaffolds to estimate their depth
conda activate samtools_v1.21

# genome.fa is the assembled scaffolds from above (myloasm or metaflye)
minimap2 -a genome.fa -o genome.sam -t 20 -x map-pb bac_reads.fastq.gz
samtools view -h -b -S genome.sam -o genome.bam --threads 20
samtools sort genome.bam -o genome.sorted.bam --threads 20
samtools index genome.sorted.bam

# 2 - binning
conda activate semibin_v2.2.0
SemiBin2 single_easy_bin \
        --environment global \
        --sequencing-type long_read \
	--processes 20 \
        -i genome.fa \
        -b genome.sorted.bam \
        -o step3_semibin_output
```

#### d - Estimating bacterial titer
To estimate the titer of bacteria, mosdepth (https://github.com/brentp/mosdepth) was used to calculate the sequencing depth of host and bacterial genomes with MAPQ >= 1. 
```
conda activate mosdepth_v0.3.8
mosdepth -f genome.fa -n -t 20 --mapq 1 mapq1 genome.sorted.bam
```

#### e - Completeness of metagenomic assembled genomes (MAGs)
We used CheckM2 to evaluate the completeness of MAGs.
```
conda activate checkm2_v1.1.0

checkm2 predict --threads 20 --input ./step3_semibin_output/output_bins/ --output-directory ./step5_checkm2_output/ --extension gz
```

### f - Species assignment
We used GTDBTk (https://github.com/Ecogenomics/GTDBTk) to assign MAGs.
```
conda activate gtdbtk_v2.5.2
gtdbtk classify_wf --genome_dir ./step3_semibin_output/output_bins/ --out_dir ./step6_gtdbtk_output/ --cpus 40 --extension gz
```

## 2 - Genome assembly and annotation of bacterial isolates
Lastly, to confirm the metagenomic assembled scaffolds, we sequenced the cultured bacteria. 
```
conda activate canu_v2.2
# -d <assembly-directory>, with output files named using the -p <assembly-prefix>
canu -p my_assembly -d canu_output genomeSize=5m -pacbio reads.fastq.gz

# Using the web server of GeneMarkS to annotate protein-coding genes.
https://exon.gatech.edu/genemarks.cgi

#Repetitive sequences


```


