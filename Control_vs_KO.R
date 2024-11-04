library(DiffBind)

setwd("/Users/swethayadavalli/Desktop/diffbind")
Control_vs_KO<-dba.peakset(NULL,
                               peaks="/Users/swethayadavalli/Desktop/SRR21677765_new_peaks_peaks.xls",
                               peak.caller="macs", sampID="Control-1",condition = "control", replicate=1, bamReads = "/Users/swethayadavalli/Desktop/SRR21677765/SRR21677765_sorted.bam")
Control_vs_KO<-dba.peakset(Control_vs_KO,
                               peaks="/Users/swethayadavalli/Desktop/SRR21677764_new_peaks_peaks.xls",
                               peak.caller="macs", sampID="Control-2",condition = "control", replicate=2, bamReads = "/Users/swethayadavalli/Desktop/SRR21677764/SRR21677764_sorted.bam")

Control_vs_KO<-dba.peakset(Control_vs_KO,
                               peaks="/Users/swethayadavalli/Desktop/SRR21677762_new_peaks_peaks.xls",
                               peak.caller="macs", sampID="KO-1",condition = "KO", replicate=1, bamReads = "/Users/swethayadavalli/Desktop/SRR21677762/SRR21677762_sorted.bam")
Control_vs_KO<-dba.peakset(Control_vs_KO,
                               peaks="/Users/swethayadavalli/Desktop/SRR21677763_new_peaks_peaks.xls",
                               peak.caller="macs", sampID="KO-2",condition = "KO", replicate=2, bamReads = "/Users/swethayadavalli/Desktop/SRR21677763/SRR21677763_sorted.bam")

Control_vs_KO_counts<-dba.count(Control_vs_KO, bParallel = TRUE, score=DBA_SCORE_READS)

png('/Users/swethayadavalli/Desktop/volcano2.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(Control_vs_KO,Control_vs_KO$masks$control, main = "Open chromatic region overlaps in Control replicates")
dev.off()

png('/Users/swethayadavalli/Desktop/KO.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(Control_vs_KO,Control_vs_KO$masks$KO, main = "Open chromatic region overlaps in the KO replicates")
dev.off()










png('/Users/swethayadavalli/Desktop/heatmap2.png', units="in", width=10, height=12, res=800)
plot(Control_vs_KO_counts, main="Correlation plot of studied samples")
dev.off()



Control_vs_KO_occupancy<-dba.peakset(Control_vs_KO, consensus = DBA_CONDITION,minOverlap = 0.33)
png('/Users/swethayadavalli/Desktop/occupancy_consensus_peaks_Control_vs_KO.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(Control_vs_KO_occupancy,Control_vs_KO_occupancy$masks$Consensus, main = "Overlap of consensus peaks of Control and KO")
dev.off()


#### Affinity analysis - differential accessibility analysis
counts<-dba.count(Control_vs_KO, bParallel = TRUE, score=DBA_SCORE_READS)
Control_vs_KO_counts<-dba.contrast(counts, categories=DBA_CONDITION,minMembers = 2)

Control_vs_KO_analyzed<-dba.analyze(Control_vs_KO_counts)
Control_vs_KO_DE_peaks <- dba.report(Control_vs_KO_analyzed)
Control_vs_KO_DE_peaks_DF<-as.data.frame(Control_vs_KO_DE_peaks)

write.table(Control_vs_KO_DE_peaks_DF[abs(Control_vs_KO_DE_peaks_DF$Fold)>1 & Control_vs_KO_DE_peaks_DF$FDR<0.05,], "/Users/swethayadavalli/Desktop/diff_regions_KO_vs_Control.txt", sep = "\t", quote = F)


library(ggplot2)
dba.plotHeatmap(Control_vs_KO_analyzed, contrast=1, bScale="row", col=brewer.pal(9, "Blues"))
dba.plotPCA(Control_vs_KO_counts)
dba.show(Control_vs_KO_analyzed)
dba.show(Control_vs_KO_DE_peaks_DF)
dba.plotMA(Control_vs_KO_analyzed)
dba.plotPCA(Control_vs_KO_analyzed)
plot(counts)


pvals <- dba.plotBox(Control_vs_KO_analyzed)
