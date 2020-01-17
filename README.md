# bmalayi_directrnaseq
 
# Set software and directory paths

## Software

```{bash, eval = F}
MINIMAP2_BIN_DIR=/usr/local/packages/minimap2-2.9/bin
SALMON_BIN_DIR=/local/aberdeen2rw/julie/Matt_dir/packages/salmon_v1.1.0/bin
SAMTOOLS_BIN_DIR=/usr/local/packages/samtools-1.9/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
```
## Directories

```{bash, eval = F}
READS_DIR=/local/aberdeen2rw/julie/Matt_dir/bmalayi_directrna
REFERENCES_DIR=/local/projects-t3/EBMAL/mchung_dir/bmalayi_directrna/references
SCRIPTS_DIR=/home/mattchung/scripts/
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/bmalayi_directrna/

MINION_FASTQ_PATH=/local/aberdeen2rw/julie/Matt_dir/native_rnaseq/fastq/fastq_runid_99e99f03c34e17a40f27574a4cccea4f0f08e63c.fastq
```

## Create directories

# Determine whether direct RNA sequencing with the MinION is biased towards short transcripts

## Set up B. malayi + wBm combined reference files

### Download B. malayi and wBm reference files from WormBase and SRA

```{bash, eval = F}
wget -O "$REFERENCES_DIR"/wBm.fna.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.fna.gz
wget -O "$REFERENCES_DIR"/bmalayi.fna.gz ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WS275.genomic.fa.gz
gunzip "$REFERENCES_DIR"/wBm.fna.gz
gunzip "$REFERENCES_DIR"/bmalayi.fna.gz

wget -O "$REFERENCES_DIR"/wBm.gff.gz ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/385/GCF_000008385.1_ASM838v1/GCF_000008385.1_ASM838v1_genomic.gff.gz
wget -O "$REFERENCES_DIR"/bmalayi.gff.gz ftp://ftp.wormbase.org/pub/wormbase/releases/WS275/species/b_malayi/PRJNA10729/b_malayi.PRJNA10729.WS275.annotations.gff3.gz
gunzip "$REFERENCES_DIR"/wBm.gff.gz
gunzip "$REFERENCES_DIR"/bmalayi.gff.gz
```

### Create nucleotide gene FASTA files

#### Input Sets
```{bash, eval = F}
## wBm
FNA="$REFERENCES_DIR"/wBm.fna
GFF3="$REFERENCES_DIR"/bmalayi.fna

## B. malayi
FNA="$REFERENCES_DIR"/wBm.gff
GFF3="$REFERENCES_DIR"/bmalayi.gff
```

#### Commands
```{bash, eval = F}
"$SCRIPTS_DIR"/createnuccdsfasta.sh -n "$FNA" -g "$GFF3" -f gene -i ID > "$(echo "$FNA" | sed "s/[.]fna$/.gene.fna/g")" &
```

### Create a combined reference from the B. malayi and wBm nucleotide gene FASTA files

#### Commands
```{bash, eval = F}
cat "$REFERENCES_DIR"/b_malayi.gene.fna "$REFERENCES_DIR"/wBm.gene.fna > "$REFERENCES_DIR"/bmalayi_wBm_combined.gene.fna
```

## Download Illumina reads from adult female B. malayi from the SRA

#### Inputs
```{bash, eval = F}
SRR_ID=SRR3111490
```

#### Commands
```{bash, eval = F}
qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" -O "$READS_DIR" --gzip
```

## Quantify B. malayi genes from Illumina sequenced adult female samples

#### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wBm.gene.fna
THREADS=4
SRR_ID=SRR3111490
```

#### Commands
```{bash, eval = F}
mkdir "$WORKING_DIR"/quant

"$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$NUC_GENE_FNA" -i "$NUC_GENE_FNA".salmon.index

echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$READS_DIR"/"$SRR_ID"_1.fastq.gz -2 "$READS_DIR"/"$SRR_ID"_2.fastq.gz -p "$THREADS" -o "$WORKING_DIR"/quant/"$SRR_ID".salmon_optimized.counts --validateMappings  --allowDovetail" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/
```

## Quantify B. malayi genes from MinION sequenced adult female samples

#### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wBm.gene.fna
THREADS=4
```

#### Commands
```{bash, eval = F}

"$MINIMAP2_BIN_DIR"/minimap2 -ax splice -uf -k14 "$NUC_GENE_FNA" "$MINION_FASTQ_PATH" |  "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$WORKING_DIR"/bmalayi_minion.bam -


echo -e ""$SALMON_BIN_DIR"/salmon quant -t "$NUC_GENE_FNA" --libType A -a "$WORKING_DIR"/bmalayi_minion.bam -p "$THREADS" -o "$WORKING_DIR"/quant/MinION.salmon_optimized.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/
```