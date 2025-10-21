# Characterization of the dynamic microbiome evolution across thrips species
A pipeline for analyzing the microbiome of thrips

## characterization of thrips microbiome

### short reads
#FastQ data were preprocessed and quality-controlled using fastp, an ultrafast all-in-one tool. https://github.com/OpenGene/fastp
#  download specified version, i.e. fastp v0.23.4
wget http://opengene.org/fastp/fastp.0.23.4
mv fastp.0.23.4 fastp
chmod a+x ./fastp
fastp -i sample1_1.fastq.gz -I sample1_2.fastq.gz -o sample1_1.clean.fastq.gz -O sample1_2.clean.fastq.gz
fastp -i sample2_1.fastq.gz -I sample2_2.fastq.gz -o sample2_1.clean.fastq.gz -O sample2_2.clean.fastq.gz
fastp -i sample3_1.fastq.gz -I sample3_2.fastq.gz -o sample3_1.clean.fastq.gz -O sample3_2.clean.fastq.gz
# Kraken2
mamba activate kraken2
time kraken2  --threads 52 --db nt --paired sample1_1.clean.fastq.gz sample1_2.clean.fastq.gz --report sample1.kreport >sample1.kraken
time kraken2  --threads 52 --db nt --paired sample2_1.clean.fastq.gz sample2_2.clean.fastq.gz --report sample2.kreport >sample2.kraken
time kraken2  --threads 52 --db nt --paired sample3_1.clean.fastq.gz sample3_2.clean.fastq.gz --report sample3.kreport >sample3.kraken

# KrakenTools is a suite of scripts to be used alongside the Kraken, KrakenUniq, Kraken 2, or Bracken programs. These scripts are designed to help Kraken users with downstream analysis of Kraken results. https://github.com/jenniferlu717/KrakenTools/

# kreport2mpa.py This program takes a Kraken report file and prints out a mpa (MetaPhlAn) -style TEXT file
python kreport2mpa.py -r sample1.kreport -o sample1.MPA.txt
python kreport2mpa.py -r sample2.kreport -o sample2.MPA.txt
python kreport2mpa.py -r sample3.kreport -o sample3.MPA.txt

# combine_mpa.py This program combines multiple outputs from kreport2mpa.py. Files to be combined must have been generated using the same kreport2mpa.py options.
python combine_mpa.py -i sample1.MPA.txt sample2.MPA.txt sample3.MPA.txt  -o 3samples.txt

# kraken2alpha.R https://github.com/YongxinLiu/EasyMicrobiome/blob/master/script/
Rscript kraken2alpha.R -i 3samples.txt

# Microbiome Helper An assortment of scripts to help process and automate various microbiome and metagenomic bioinformatic tools.https://github.com/LangilleLab/microbiome_helper/tree/master 
perl metaphlan_to_stamp.pl 3samples.norm.txt >3samples.norm.txt.spf
Rscript metaphlan_hclust_heatmap.R -i 3samples.norm.txt.spf -t Species -n 25 -o heatmap_Species


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
