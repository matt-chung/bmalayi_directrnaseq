# bmalayi_directrnaseq
 
# Set software and directory paths

## Software

```{bash, eval = F}
SALMON_BIN_DIR=/usr/local/packages/salmon-0.13.1/bin
SRATOOLKIT_BIN_DIR=/usr/local/packages/sratoolkit-2.9.0/bin
```
## Directories

```{bash, eval = F}
READS_DIR=/local/aberdeen2rw/julie/Matt_dir/bmalayi_directrna
REFERENCES_DIR=/local/projects-t3/EBMAL/mchung_dir/bmalayi_directrna/references
SCRIPTS_DIR=/home/mattchung/scripts/
WORKING_DIR=/local/projects-t3/EBMAL/mchung_dir/bmalayi_directrna/
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

#### Input Sets:
```{bash, eval = F}
## wBm
FNA="$REFERENCES_DIR"/wBm.fna
GFF3="$REFERENCES_DIR"/bmalayi.fna

## B. malayi
FNA="$REFERENCES_DIR"/wBm.gff
GFF3="$REFERENCES_DIR"/bmalayi.gff
```

#### Commands:
```{bash, eval = F}
"$SCRIPTS_DIR"/createnuccdsfasta.sh -n "$FNA" -g "$GFF3" -f gene -i ID > "$(echo "$FNA" | sed "s/[.]fna$/.gene.fna/g")" &
```

### Create a combined reference from the B. malayi and wBm nucleotide gene FASTA files

#### Commands:
```{bash, eval = F}
cat "$REFERENCES_DIR"/b_malayi.gene.fna "$REFERENCES_DIR"/wBm.gene.fna > "$REFERENCES_DIR"/bmalayi_wBm_combined.gene.fna
```

## Download Illumina reads from adult female B. malayi from the SRA

```{bash, eval = F}
SRR_ID=SRR3111490
```

```{bash, eval = F}
qsub -P jdhotopp-lab -l mem_free=2G -N fastq_dump -wd "$READS_DIR" -b y "$SRATOOLKIT_BIN_DIR"/fastq-dump --split-files "$SRR_ID" -O "$READS_DIR"
```

## Download Illumina reads from adult female B. malayi from the SRA
