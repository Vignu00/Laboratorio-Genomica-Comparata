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
blobtools view -i <JSON_FILE> -o <OUTPUT_PREFIX>
blobtools plot -i <JSON_FILE> -o <OUTPUT_PREFIX>
blobtools view -i Anoste_blob.blobDB.json -o Anoste
blobtools plot -i Anoste_blob.blobDB.json -o Anoste
```
grep -v '^#' Anoste.Anoste_blob.blobDB.table.txt | awk '$6 != "Arthropoda" {print $1}' > contaminants_contigs.txt #Selezionare tutte le contigs contaminanti (non artropoda) e raccoglierle in un file che servirà come file patterns.txt nel comando successivo awk
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' in.fa | grep -w -v -Ff patterns.txt - | tr "\t" "\n" > out.fa
awk '{ if ((NR>1)&&($0~/^>/)) { printf("\n%s", $0); } else if (NR==1) { printf("%s", $0); } else { printf("\t%s", $0); } }' Anoste_pol.fasta | grep -w -v -Ff contaminants_contigs.txt - | tr "\t" "\n" > Anoste_pol_nocont.fa # Creato file fasta senza contaminanti
cd ..
mkdir 03_scaffolding (/Laboratorio-Genomica-Comparata/03_GenomeAssembly/03_scaffolding)
ln -s ../02_contaminants/Anoste_pol_nocont.fa
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/01_Anoste_reference/GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna
conda activate assembly # Qui c'è ragtag

# Reference-based scaffolding

```bash
ragtag.py correct -t 20 <REFERENCE_GENOME> <DRAFT_GENOME>
ragtag.py scaffold -C -t 20 -o <OUTPUT_DIR> <REFERENCE_GENOME> <CORRECTED_DRAFTGENOME>
ragtag.py correct -t 10 GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna Anoste_pol_nocont.fa
ragtag.py scaffold -C -t 10 -o Anoste_chr GCF_013141755.1_UCI_ANSTEP_V1.0_genomic.fna ragtag_output/ragtag.correct.fasta
```
cd Anoste_chr/
grep ">" ragtag.scaffold.fasta | head # Selezionare le prime 3 righe (veri cromosomi)
cd ..
nano chromosome.txt # File di patterns dentro incollare le 3 righe >NC_050201.1_RagTag >NC_050202.1_RagTag >NC_050203.1_RagTag
grep -A1 -Ff chromosome.txt Anoste_chr/ragtag.scaffold.fasta > Anoste_chr.fasta

# Genome annotation

mkdir 04_GenomeAnnotation
cd 04_GenomeAnnotation/
mkdir 00_round1
mkdir 01_round2
mkdir 02_round3
cd 00_round1/
ln -s ../../03_GenomeAssembly/03_scaffolding/Anoste_chr.fasta
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation/GCF_013141755.1_UCI_ANSTEP_V1.0_protein.faa Anoste.faa
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation/GCF_943734735.2_idAnoGambNW_F1_1_protein.faa Anogam.faa
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation/GCF_000001215.4_Release_6_plus_ISO1_MT_protein.faa Dromel.faa
ln -s /home/PERSONALE/mirko.martini3/Lab_CompGeno/00_practice/00_data/02_Annotation/Anoste_RepeatModeler_library.fa
conda activate assembly
maker -CTL
nano maker_opts.ctl # Dentro modificare: 

# maker_opts.ctl
```bash
#-----Genome (these are always required)
genome=Anoste_chr.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff= #aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein=Anogam.faa,Anoste.faa,Dromel.faa #protein sequence file in fasta format (i.e. from mutiple organisms)
protein_gff=  #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for DFam masking in RepeatMasker
rmlib=Anoste_RepeatModeler_library.fa #provide an organism specific repeat library in fasta format for RepeatMasker #Fatto library da Mirko
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
run_evm=0 #run EvidenceModeler, 1 = yes, 0 = no
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
snoscan_meth= #-O-methylation site fileto have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no
allow_overlap= #allowed gene overlap fraction (value from 0 to 1, blank for default)

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=9 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=1 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=1 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=50 #require at least this many amino acids in predicted proteins
alt_splice=1 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=0 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
min_intron=20 #minimum intron length (used for alignment polishing)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
```



