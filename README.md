# Rep_Labgenomica

# kmer based survey

```bash
# Comando per fare il QC delle raw reads appena sequenziate
fastqc /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672503_1.fastq.gz /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672503_2.fastq.gz
```
```bash
# Comando per trimming delle sequenze a cui è stato fatto il QC
trimmomatic PE -threads 8 -phred33 SRR11672503_1.fastq.gz SRR11672503_2.fastq.gz SRR11672503_1_paired.fastq SRR11672503_1_unpaired.fastq SRR11672503_2_paired.fastq SRR11672503_2_unpaired.fastq ILLUMINACLIP:/opt/miniforge3/envs/assembly/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> stats_trimmomatic.log
```
```bash
conda activate kat
# Comando per analizzare kmers
kat hist -t <THREADS> -m <MER LENGTH> -o <OUTPUT PREFIX> <paired_1> <paired_2>
kat hist -t 8 -m 27 -o anoste_27 SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq  
```

# Assembly

tmux a -s genome_assembly
mkdir 03_GenomeAssembly
cd 03_GenomeAssembly

```bash
# Comando per layout
wtdbg2 -x rs -g <EXPECTED_GENOMESIZE> -t <NUMBERS_OF_CORES> -i <FASTQ> -fo <OUT_PREFIX> #Default parameters
wtdbg2 -x rs -g <EXPECTED_GENOMESIZE> -t 9 -i SRR11672506.fastq.gz -fo Anoste_raw.fasta
```

```bash
# Comando per consensus
wtpoa-cns -t <NUMBER_OF_CORES> -i <LAYOUT_FILE-.ctg.lay.gz> -fo <OUT_PREFIX>
wtpoa-cns -t 9 -i Anoste_raw.ctg.lay.gz> -fo Anoste_raw.fasta
```
mkdir 00_Assembly_raw
cd 00_Assembly_raw
mv * 00_Assembly_raw
coda activate sequence #Enviroment in cui c'è BUSCO
export NUMEXPR_MAX_THREADS=80

```bash
# Comando per BUSCO
busco -m <MODE> -l <LINEAGE> -c <CPU_NUMBER> -o <OUTPUT_NAME> -i <INPUT>
busco -m geno -l $BUSCO/diptera_odb12 -c 8 -o Anoste_raw -i Anoste_raw.fasta
```
tmux new -s mapping
cd Laboratorio-Genomica-Comparata/03_GenomeAssembly/
mkdir 01_polishing
cd 01_polishing/
nano mapping.sh # Creare e aprire file e incollare comandi (script)
#!/bin/bash # Dentro file creato

# Mapping short and long reads on assembly

"$CPU"=$1
"$INPUT"=$2
"$SR_1"=$3
"$SR_2"=$4
"$PREFIX"=$5

# Mapping short and long reads on assembly
minimap2 -ax <PRESET_Options> --MD -t <cores_number> <target.fa> <fastq1> <fastq2> > <SAM_FILE>
minimap2 -ax sr --MD -t "$CPU" "$INPUT" "$SR_1" "$SR_2" > "$PREFIX".sam

# Convertion SAM into BAM
samtools view -b <SAM FILE> > <BAM FILE>
samtools view -b "$PREFIX".sam > "$PREFIX".bam

# Rimuovere file SAM
rm <SAM FILE>
rm "$PREFIX".sam

# Sorting and indexing of BAM file
samtools sort -@<NUM_THREADS> -o <OUT_FILE> <INFILE>
samtools index <INFILE>
samtools sort -@"$CPU" -o "$PREFIX"_sorted.bam "$PREFIX".bam
samtools index "$PREFIX"_sorted.bam

ln -s ../00_Assembly_raw/Anoste_raw.fasta .
ln -s ../../02_kmersurvey/SRR11672503_1_paired.fastq .
ln -s ../../02_kmersurvey/SRR11672503_2_paired.fastq .
conda activate assembly # SAM tool
bash mapping.sh 8 Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq Anoste_pol



