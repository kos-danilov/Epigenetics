# 19. Positional anotation with ChIPseeker
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
peaks <- readPeakFile("~/chipseq/span/14_hg19_0.05_5.peak")
peaks <- readPeakFile("~/chipseq/sicer/14_hg19-W200-G600-FDR0.05-island.bed")

peakAnno <- annotatePeak(peaks, tssRegion=c(-3000, -3000), TxDb=txdb, annoDb="org.Hs.eg.db")
plotAnnoPie(peakAnno)
plotAnnoBar(peakAnno)
plotDistToTSS(peakAnno, title="Distribution wrt to TSS")

# 20. Peaks annotation
library(ChIPpeakAnno)
# Correct format
gr1 <- toGRanges("~/chipseq/span/14_hg19_0.05_5.peak3", format="BED", header=FALSE)

# 21. Positional annotation of replicated peaks
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
genomicElementDistribution(gr1,
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))

# 22. Genomic annotation as UpsetR plot
library(UpSetR)
x <- genomicElementUpSetR(gr1,
                          TxDb.Hsapiens.UCSC.hg19.knownGene)
upset(x$plotData, nsets=13, nintersects=NA)

# 23. Gene set enrichment of peaks
gr1.anno <- annotatePeakInBatch(gr1,
                                     AnnotationData=annoData,
                                     output="nearestBiDirectionalPromoters",
                                     bindingRegion=c(-2000, 500))
head(gr1.anno)

# 24. Annotate peaks with entrez genes, and run GO enrichment
library(org.Hs.eg.db)
gr1.anno <- addGeneIDs(gr1.anno, "org.Hs.eg.db", IDs2Add = "entrez_id")
head(overlaps.anno)

library("DBI")
over <- getEnrichedGO(gr1.anno, orgAnn="org.Hs.eg.db", condense=TRUE)
enrichmentPlot(over)

# 26. Enrichment vs Reactome DB
library(reactome.db)
path <- getEnrichedPATH(gr1.anno, "org.Hs.eg.db", "reactome.db", maxP=.05)
enrichmentPlot(path)

