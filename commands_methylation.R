# 20. Regions overlap
library(ChIPpeakAnno)
# Correct format
donor1_1 <- toGRanges("/home/student/bs_seq/meth_CpG_hypermr/DC552_HI48.merged_CpG.hypermr", format="BED", header=FALSE)
donor1_2 <- toGRanges("/home/student/bs_seq/meth_CpG_hypermr/DC552_NI48.merged_CpG.hypermr", format="BED", header=FALSE)

donor2_1 <- toGRanges("/home/student/bs_seq/meth_CpG_hypermr/DC555_HI48.merged_CpG.hypermr", format="BED", header=FALSE)
donor2_2 <- toGRanges("/home/student/bs_seq/meth_CpG_hypermr/DC555_NI48.merged_CpG.hypermr", format="BED", header=FALSE)

ol <- findOverlapsOfPeaks(donor1_1, donor1_2, donor2_1, donor2_2)

makeVennDiagram(ol)

donor1_1 <- toGRanges("/home/student/bs_seq/meth_CpG_hmr/DC552_HI48.merged_CpG.hmr", format="BED", header=FALSE)
donor1_2 <- toGRanges("/home/student/bs_seq/meth_CpG_hmr/DC552_NI48.merged_CpG.hmr", format="BED", header=FALSE)

donor2_1 <- toGRanges("/home/student/bs_seq/meth_CpG_hmr/DC555_HI48.merged_CpG.hmr", format="BED", header=FALSE)
donor2_2 <- toGRanges("/home/student/bs_seq/meth_CpG_hmr/DC555_NI48.merged_CpG.hmr", format="BED", header=FALSE)

ol <- findOverlapsOfPeaks(donor1_1, donor1_2, donor2_1, donor2_2)

makeVennDiagram(ol)

# Positional annotation
# can not install TxDb.Hsapiens.UCSC.hg38.knownGene on server
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
peaks <- GRangesList(rep1=donor1_1, rep2=donor1_2, rep3=donor2_1, rep4=donor2_2)
genomicElementDistribution(peaks,
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))

#Distribution vs genomic features
up=toGRanges("/home/student/bs_seq/meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3.up.bed3", format="BED", header=FALSE)
down=toGRanges("/home/student/bs_seq/meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3.down.bed3", format="BED", header=FALSE)
combo=toGRanges("/home/student/bs_seq/meth_CpG_dmrs/dmrs_avgdiff_0.025_ncyto_3.bed3", format="BED", header=FALSE)
genomicElementDistribution(up,
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))
genomicElementDistribution(down,
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))
genomicElementDistribution(combo,
                           TxDb = TxDb.Hsapiens.UCSC.hg19.knownGene,
                           promoterRegion=c(upstream=2000, downstream=500),
                           geneDownstream=c(upstream=0, downstream=2000))


