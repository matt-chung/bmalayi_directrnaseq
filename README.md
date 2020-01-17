# bmalayi_directrnaseq
 
<!-- MarkdownTOC -->

- [Set software and directory paths](#set-software-and-directory-paths)
  - [Software](#software)
  - [Directories](#directories)
  - [Create directories](#create-directories)
- [Determine whether direct RNA sequencing with the MinION is biased towards short transcripts](#determine-whether-direct-rna-sequencing-with-the-minion-is-biased-towards-short-transcripts)
  - [Set up B. malayi + wBm combined reference files](#set-up-b-malayi--wbm-combined-reference-files)
    - [Download B. malayi and wBm reference files from WormBase and SRA](#download-b-malayi-and-wbm-reference-files-from-wormbase-and-sra)
    - [Create nucleotide gene FASTA files](#create-nucleotide-gene-fasta-files)
      - [Input Sets](#input-sets)
      - [Commands](#commands)
    - [Create a combined reference from the B. malayi and wBm nucleotide gene FASTA files](#create-a-combined-reference-from-the-b-malayi-and-wbm-nucleotide-gene-fasta-files)
      - [Commands](#commands-1)
  - [Download Illumina reads from adult female B. malayi from the SRA](#download-illumina-reads-from-adult-female-b-malayi-from-the-sra)
      - [Inputs](#inputs)
      - [Commands](#commands-2)
  - [Quantify B. malayi genes from Illumina sequenced adult female samples](#quantify-b-malayi-genes-from-illumina-sequenced-adult-female-samples)
    - [Create index for Salmon read-mode quantification](#create-index-for-salmon-read-mode-quantification)
      - [Inputs](#inputs-1)
      - [Commands](#commands-3)
    - [Quantify genes from Illumina reads using read-mode Salmon](#quantify-genes-from-illumina-reads-using-read-mode-salmon)
      - [Inputs](#inputs-2)
      - [Commands](#commands-4)
  - [Quantify B. malayi genes from MinION sequenced adult female samples](#quantify-b-malayi-genes-from-minion-sequenced-adult-female-samples)
    - [Map MinION reads to gene fasta for alignment-mode Salmon quantification](#map-minion-reads-to-gene-fasta-for-alignment-mode-salmon-quantification)
      - [Inputs](#inputs-3)
      - [Commands](#commands-5)
    - [Quantify genes from MinION reads using alignment-mode Salmon](#quantify-genes-from-minion-reads-using-alignment-mode-salmon)
      - [Inputs](#inputs-4)
      - [Commands](#commands-6)
  - [Compare the ratio of normalized counts from the MinION and Illumina sequenced samples](#compare-the-ratio-of-normalized-counts-from-the-minion-and-illumina-sequenced-samples)
    - [Set R inputs](#set-r-inputs)
    - [Load R packages and view sessionInfo](#load-r-packages-and-view-sessioninfo)
    - [Load R functions](#load-r-functions)
    - [Create normalized counts data frame](#create-normalized-counts-data-frame)
    - [Exclude Wolbachia genes from counts dataframe](#exclude-wolbachia-genes-from-counts-dataframe)
    - [Calculates log2 MinION:Illumina normalized count ratio for each gene](#calculates-log2-minionillumina-normalized-count-ratio-for-each-gene)
    - [Identifies the number of genes with a log2 MinION:Illumina normalized count ratio >2 and <-2](#identifies-the-number-of-genes-with-a-log2-minionillumina-normalized-count-ratio-2-and--2)
    - [Plot gene length versus log2 MinION:Illumina normalized count ratios for each gene](#plot-gene-length-versus-log2-minionillumina-normalized-count-ratios-for-each-gene)
    - [Use density plot to visualize the gene length versus log2 MinION:Illumina normalized count ratios for each gene](#use-density-plot-to-visualize-the-gene-length-versus-log2-minionillumina-normalized-count-ratios-for-each-gene)

<!-- /MarkdownTOC -->


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
OUTPUT_DIR=/local/projects-t3/EBMAL/mchung_dir/bmalayi_directrna/output

MINION_FASTQ_PATH=/local/aberdeen2rw/julie/Matt_dir/native_rnaseq/fastq/fastq_runid_99e99f03c34e17a40f27574a4cccea4f0f08e63c.fastq
```

## Create directories
```{bash, eval = F}
mkdir "$WORKING_DIR"/quant
```

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

### Create index for Salmon read-mode quantification

#### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wBm.gene.fna
```

#### Commands
```{bash, eval = F}
"$SALMON_BIN_DIR"/salmon index --keepDuplicates -t "$NUC_GENE_FNA" -i "$NUC_GENE_FNA".salmon.index
```

### Quantify genes from Illumina reads using read-mode Salmon

#### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wBm.gene.fna
THREADS=4
SRR_ID=SRR3111490
```

#### Commands
```{bash, eval = F}
echo -e ""$SALMON_BIN_DIR"/salmon quant -i "$NUC_GENE_FNA".salmon.index --libType A -1 "$READS_DIR"/"$SRR_ID"_1.fastq.gz -2 "$READS_DIR"/"$SRR_ID"_2.fastq.gz -p "$THREADS" -o "$WORKING_DIR"/quant/"$SRR_ID".salmon_optimized.counts --validateMappings  --allowDovetail" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/
```

## Quantify B. malayi genes from MinION sequenced adult female samples

### Map MinION reads to gene fasta for alignment-mode Salmon quantification

Read-mode Salmon quantification was run using the MinION reads but returned 0 reads mapped.

#### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wBm.gene.fna
THREADS=4
```

#### Commands
```{bash, eval = F}
"$MINIMAP2_BIN_DIR"/minimap2 -ax splice -uf -k14 "$NUC_GENE_FNA" "$MINION_FASTQ_PATH" |  "$SAMTOOLS_BIN_DIR"/samtools view -bhSo "$WORKING_DIR"/bmalayi_minion.bam -
```

### Quantify genes from MinION reads using alignment-mode Salmon

#### Inputs
```{bash, eval = F}
NUC_GENE_FNA="$REFERENCES_DIR"/bmalayi_wBm.gene.fna
THREADS=4
```

#### Commands
```{bash, eval = F}
echo -e ""$SALMON_BIN_DIR"/salmon quant -t "$NUC_GENE_FNA" --libType A -a "$WORKING_DIR"/bmalayi_minion.bam -p "$THREADS" -o "$WORKING_DIR"/quant/MinION.salmon_optimized.counts" | qsub -P jdhotopp-lab -q threaded.q -pe thread "$THREADS" -l mem_free=5G -N salmon -wd "$WORKING_DIR"/quant/
```

## Compare the ratio of normalized counts from the MinION and Illumina sequenced samples

### Set R inputs

```{r}
QUANT_DIR="Z:/EBMAL/mchung_dir/bmalayi_directrna/quant"
```

### Load R packages and view sessionInfo

```{r}
library(ggplot2)
library(reshape2)

sessionInfo()
```

```{r, eval = F}
R version 3.5.1 (2018-07-02)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 7 x64 (build 7601) Service Pack 1

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] see_0.3.0          rlist_0.4.6.1      reshape2_1.4.3     pvclust_2.0-0      matrixStats_0.55.0 ggplot2_3.2.1     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.1        pillar_1.3.1      compiler_3.5.1    plyr_1.8.4        parameters_0.3.0  tools_3.5.1       tibble_2.1.1     
 [8] gtable_0.3.0      pkgconfig_2.0.2   rlang_0.4.0       rstudioapi_0.10   yaml_2.2.0        xfun_0.6          gridExtra_2.3    
[15] withr_2.1.2       bayestestR_0.4.0  dplyr_0.8.3       stringr_1.4.0     knitr_1.22        grid_3.5.1        tidyselect_0.2.5 
[22] glue_1.3.1        data.table_1.12.2 R6_2.4.0          purrr_0.3.2       magrittr_1.5      ggridges_0.5.1    scales_1.0.0     
[29] assertthat_0.2.1  insight_0.7.1     effectsize_0.0.1  colorspace_1.4-1  labeling_0.3      stringi_1.4.3     lazyeval_0.2.2   
[36] munsell_0.5.0     crayon_1.3.
```

### Load R functions
```{r}
lm_eqn <- function(df){
    x <- df[,1]
    y <- df[,2]
    m <- lm(y ~ x, df);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}
```

### Create normalized counts data frame
```{r}
illumina_quant <- read.delim(paste0(QUANT_DIR,"/SRR3111490.salmon_optimized.counts/quant.sf"))
minion_quant <- read.delim(paste0(QUANT_DIR,"/MinION.salmon_optimized.counts/quant.sf"))

minion_quant <- minion_quant[match(rownames(illumina_quant),rownames(minion_quant)),]

counts <- as.data.frame(cbind(illumina_quant[,5]/sum(illumina_quant[,5]),
                              minion_quant[,5]/sum(minion_quant[,5])))
colnames(counts) <- c("Illumina","MinION")
rownames(counts) <- illumina_quant[,1]
```

### Exclude Wolbachia genes from counts dataframe
```{r}
counts <- counts[!grepl("^gene",rownames(counts)),]
```

### Calculates log2 MinION:Illumina normalized count ratio for each gene

Genes with a ratio of >2 have significantly more counts from MinION while genes with a ratio of <-2 have significantly more counts from Illumina.

```{r,fig.height=4,fig.width=6}
counts.plot.df <- as.data.frame(cbind(counts$MinION/counts$Illumina,
                                      illumina_quant[!grepl("^gene",illumina_quant[,1]),2]))
counts.plot.df[,1] <- log2(counts.plot.df[,1])
```

### Identifies the number of genes with a log2 MinION:Illumina normalized count ratio >2 and <-2

The number of genes with 4x greater counts in the MinION sample:
```{r}
table(counts.plot.df[,1] > 2)
```

```{r, eval = F}
FALSE  TRUE 
10113  1219
```

The number of genes with 4x greater counts in the Illumina sample:
```{r}
table(counts.plot.df[,1] < -2)
```

```{r, eval = F}
FALSE  TRUE 
 7053  4279 
```

### Plot gene length versus log2 MinION:Illumina normalized count ratios for each gene

The only genes skewed towards the MinION-sequenced sample are those of a smaller gene length, indicating most of the counts from the MinION sample are being assigned to smaller genes. On the other hand, more of the larger genes are skewed towards the Illumina sample, indicating the counts are more evenly distributed among the different genes.

```{r,fig.height=5,fig.width=6}
counts.plot.df <- counts.plot.df[counts.plot.df != Inf,]
counts.plot.df <- counts.plot.df[counts.plot.df != -Inf,]
counts.plot.df <- counts.plot.df[!is.na(counts.plot.df[,1]),]

ratio.plot <- ggplot(mapping=aes(x=counts.plot.df[,2],y=counts.plot.df[,1]))+
  geom_hline(mapping=aes(yintercept=2),color="red",lty="dashed")+
  geom_hline(mapping=aes(yintercept=-2),color="red",lty="dashed")+
  geom_point(size=0.5)+
  labs(x="gene length (bp)", y = "log2 normalized MinION counts:normalized Illumina counts")+
  coord_cartesian(ylim=c(-11,11))+
  scale_x_continuous(lim=c(0,NA))+
  geom_smooth(method='lm',formula=y~x)+
  geom_text(mapping=aes(x=0, y = -10), label = lm_eqn(counts.plot.df),hjust=0,parse = TRUE)+
  theme_bw()

pdf(paste0(OUTPUT_DIR,"/ratio_plot.pdf"),
    height=5,
    width=6)
print(ratio.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/ratio_plot.png"),
    height=5,
    width=6,
    units = "in",res=300)
print(ratio.plot)
dev.off()

print(ratio.plot)
```

![Image description](/images/ratio_plot.png)

### Use density plot to visualize the gene length versus log2 MinION:Illumina normalized count ratios for each gene

```{r,fig.height=5,fig.width=6}
density_ratio.plot <- ggplot(mapping=aes(x=counts.plot.df[,2],y=counts.plot.df[,1]))+
  geom_hline(mapping=aes(yintercept=2),color="red",lty="dashed")+
  geom_hline(mapping=aes(yintercept=-2),color="red",lty="dashed")+
  stat_density_2d(aes(fill = stat(nlevel)), geom = "polygon")+
  labs(x="gene length (bp)", y = "log2 normalized MinION counts:normalized Illumina counts")+
  coord_cartesian(ylim=c(-11,11))+
  scale_x_continuous(lim=c(0,NA))+
  geom_smooth(method='lm',formula=y~x)+
  geom_text(mapping=aes(x=0, y = -10), label = lm_eqn(counts.plot.df),hjust=0,parse = TRUE)+
  scale_fill_viridis_c()+
  guides(fill = F)+
  theme_bw()

pdf(paste0(OUTPUT_DIR,"/density_ratio_plot.pdf"),
    height=5,
    width=6)
print(density_ratio.plot)
dev.off()

png(paste0(OUTPUT_DIR,"/density_ratio_plot.png"),
    height=5,
    width=6,
    units = "in",res=300)
print(density_ratio.plot)
dev.off()

print(density_ratio.plot)
```

![Image description](/images/density_ratio_plot.png)