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

tmux a -t genome_assembly
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
tmux new -t mapping
cd Laboratorio-Genomica-Comparata/03_GenomeAssembly/
mkdir 01_polishing
cd 01_polishing/
nano mapping.sh # Creare e aprire file e incollare comandi (script)


# Mapping short and long reads on assembly

```bash
#!/bin/bash # Dentro file creato
CPU=$1
INPUT=$2
SR_1=$3
SR_2=$4
PREFIX=$5
LR=$6

# Mapping SHORT reads
# Mapping short on assembly
minimap2 -ax sr --MD -t "$CPU" "$INPUT" "$SR_1" "$SR_2" > "$PREFIX"_sr.sam

# Convertion SAM into BAM
samtools view -b "$PREFIX"_sr.sam > "$PREFIX"_sr.bam

# Rimuovere file SAM
rm "$PREFIX"_sr.sam

# Sorting and indexing of BAM file
samtools sort -@"$CPU" -o "$PREFIX"_sr_sorted.bam "$PREFIX"_sr.bam
samtools index "$PREFIX"_sr_sorted.bam

# Mapping LONG reads
# Mapping long on assembly
minimap2 -ax map-pb --MD -t "$CPU" "$INPUT" "$LR" > "$PREFIX"_lr.sam

# Convertion SAM into BAM
samtools view -b "$PREFIX"_lr.sam > "$PREFIX"_lr.bam

# Rimuovere file SAM
rm "$PREFIX"_lr.sam

# Sorting and indexing of BAM file
samtools sort -@"$CPU" -o "$PREFIX"_lr_sorted.bam "$PREFIX"_lr.bam
samtools index "$PREFIX"_lr_sorted.bam
```

ln -s ../00_Assembly_raw/Anoste_raw.fasta .
ln -s ../../02_kmersurvey/SRR11672503_1_paired.fastq . # Link short reads
ln -s ../../02_kmersurvey/SRR11672503_2_paired.fastq .
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/00_reads/SRR11672506.fastq.gz . # Link long reads
conda activate assembly # Si trova qui SAM tool

```bash
# Attivare script
bash mapping.sh 8 Anoste_raw.fasta SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq Anoste_pol SRR11672506.fastq.gz
```
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/01_polishing/Anoste_raw_pb_sorted.bam .
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/01_polishing/Anoste_raw_pb_sorted.bam.bai .
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/01_polishing/Anoste_raw_sr_sorted.bam .
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/01_polishing/Anoste_raw_sr_sorted.bam.bai .
rm Anoste_raw.fasta # Non copia link se c'è file con stesso nome
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/00_assembly_raw/Anoste_raw.fasta .

```bash
# Calcolo coverage sulle short reads
mosdepth -n --fast-mode --by 500 sr_coverage Anoste_raw_sr_sorted.bam
zcat sr_coverage.regions.bed.gz | awk '{sum += $4;count++} END {print sum / count}' > coverage_sr
```
nano polishing.sh # Creo nuovo script per polishing e controllo qualità

# Polishing

R1="/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/02_KmerBased_GenomeSurvey/SRR11672503_1_paired.fastq"
R2="/home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/02_KmerBased_GenomeSurvey/SRR11672503_2_paired.fastq"
echo -e "$R1\n$R2" > reads.txt 

```bash # Dentro file creato (Alla fine non abbiamo usato script perchè QC lanciato successivamente)
hypo -d <DRAFT_Contigs> -r @<READS_PATH> -s <APPROXIMATE_GENOMESIZE> -c <SHORT_READSCOVERAGE> -b <SORTED_BAM_SR> -B <SORTED_BAM_PB> -t <NUMBER_THREADS>
hypo -d Anoste_raw.fasta -r @reads.txt -s 227054799 -c 136 -b Anoste_raw_sr_sorted.bam -B Anoste_raw_pb_sorted.bam -t 8
```
mv hypo_Anoste_raw.fasta Anoste_pol.fasta # Cambiare denominazione del risultato di hypo

```bash
# BUSCO dopo polishing 
conda activate sequence 
export NUMEXPR_MAX_THREADS=80
fold -w 80 Anoste_pol.fasta > Anoste_pol_multi.fasta # Per convertire file fasta in multi line di 80 (BUSCO non accetta single line)
rm Anoste_pol.fasta
mv Anoste_pol_multi.fasta Anoste_pol.fasta
busco -m geno -l $BUSCO/diptera_odb12 -c 9 -o Anoste_pol -i Anoste_pol.fasta
```
```bash
# Analisi QC kat
conda activate kat
kat comp -t 9 -o Anoste_pol_kat 'SRR11672503_1_paired.fastq SRR11672503_2_paired.fastq' Anoste_pol.fasta
```
cd..
mkdir 02_contaminants
cd 02_contaminants/
ln -s ../01_polishing/Anoste_pol.fasta .
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/02_contaminants/Anoste_pol_sorted.bam . # Link del re-mapping e annotazione tassonomica i quali non lanciamo 
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/02_contaminants/Anoste_pol_sorted.bam.bai .
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/03_GenomeAssembly/02_contaminants/Anoste_blast.tsv .

# Re-mapping short-reads (completa con dati da solo) (non lanciato)
```bash
minimap2 --secondary=no --MD -ax sr -t <NUMBER_CORE> <ASSEMBLY> <FASTQ_R1> <FASTQ_R2> | samtools view -Sb - > <OUT_BAMFILE>.bam
samtools sort -@<n_cores> -o <OUT_SORTED_BAMFILE>_sorted.bam Anoste_pol.bam
rm Anoste_pol.bam
samtools index Anoste_pol_sorted.bam
```

# Taxonomic annotation of the contigs (non lanciato)
```bash
blastn -query <ASSEMBLY> -db <PATH/TO/nt/> -outfmt '6 qseqid staxids bitscore std sscinames sskingdoms stitle' -max_target_seqs 25 -max_hsps 1 -num_threads 25 -evalue 1e-25 -out <OUTFILE>
```

# Construction a BlobDB data structure and assembly visualisation
```bash
blobtools create -i Anoste_pol.fasta -b Anoste_pol_sorted.bam -t Anoste_blast.tsv -o Anoste_blob # Sequence enviroment

