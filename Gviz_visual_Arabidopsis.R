## after https://support.bioconductor.org/p/90392/
library(ggbio)
library(Gviz)
library(biomaRt)
library(BSgenome.Athaliana.TAIR.TAIR9)
## all packages are within Bioconductor
## e.g.
## source("https://bioconductor.org/biocLite.R")
## if (!requireNamespace("BiocManager", quietly = TRUE))
## install.packages("BiocManager")
## BiocManager::install("ggbio")
## BiocManager::install("Gviz")
## BiocManager::install("biomaRt")
## BiocManager::install("BSgenome.Athaliana.TAIR.TAIR9")

## for Data Track (histogram of SNP numbers at the 5 nucleotide sites)
## prepared a dataset file called ‘SNPnumbers2.csv’ with proportions of SNPs at 5 nucleotide positions

SNPscores<-read.csv("SNPnumbers2.csv",header=TRUE)

## converted dataset into a GRanges object

SNPscoresgr<-makeGRangesFromDataFrame(SNPscores,keep.extra.columns=TRUE)

## [genome annotation?] data frame for ideogram track, modified from https://support.bioconductor.org/p/90392/
idTable<-data.frame(chrom=tolower(seqnames(BSgenome.Athaliana.TAIR.TAIR9)),
                    chromStart=0,chromEnd=seqlengths(BSgenome.Athaliana.TAIR.TAIR9),
                    name=seqnames(BSgenome.Athaliana.TAIR.TAIR9),
                    gieStain="gneg",stringsAsFactors=FALSE)

## ideogram track - an idealized representation of a single chromosome
idTrack<-IdeogramTrack(bands=idTable,genome="TAIR9")
mart=useEnsembl(biomart="plants_mart",host="plants.ensembl.org",dataset="athaliana_eg_gene")
bm<-useMart(host="plants.ensembl.org",biomart="plants_mart",dataset="athaliana_eg_gene")
fm<-Gviz:::.getBMFeatureMap()
fm["cdsl"]<-"cds_start"
fm["symbol"]<-"ensembl_transcript_id"

biomTrack<-BiomartGeneRegionTrack(name="ENSEMBL",symbol="LHY",biomart=bm,featureMap=fm)

gTrack<-GenomeAxisTrack()

aTrack.groups<-AnnotationTrack(start=c(37569,37373,37062),
                               width=c(272,26,142),chromosome="chr1",strand="-",group=c("exon1","exon2","exon3"),
                               genome="TAIR9",name="5'UTR",groupAnnotation="group", just.group="above")

## to get brick and stick exon-intron structure (introns as single lines and no annotation

aTrack.stacked<-AnnotationTrack(start=c(20000,37569,37373,37062,40000),
                                width=c(1,272,26,142,1),chromosome="chr1",strand=rep(c("-","*","-"),c(1,3,1)),
                                group=rep(c("5'","exon1-exon2-exon3","3'"),c(1,3,1)),
                                genome="TAIR9",name="5'UTR",groupAnnotation="group",stacking="dense")

dTrack<-DataTrack(SNPscoresgr,groups=c("T","G","C","A","N"),col=c("brown3","royalblue1","palegreen3","grey","black"),
                  type="histogram",stackedBars=TRUE,baseline=0,col.baseline="grey",name="SNP abundance")

png("LHY_5UTR_haplotype.png")

plotTracks(list(idTrack,gTrack,biomTrack,aTrack.stacked,dTrack),transcriptAnnotation="transcript",reverseStrand=TRUE,from=37020,to=37965,
                             background.panel="#FFFEDB",background.title="darkblue")

dev.off()
