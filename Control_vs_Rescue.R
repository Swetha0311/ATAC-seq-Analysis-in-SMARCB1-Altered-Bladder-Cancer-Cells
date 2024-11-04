if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DiffBind")
library(DiffBind)

setwd("/Users/swethayadavalli/Desktop/diffbind")
Control_vs_Rescue<-dba.peakset(NULL,
                          peaks="/Users/swethayadavalli/Desktop/SRR21677765_new_peaks_peaks.xls",
                          peak.caller="macs", sampID="Control-1",condition = "control", replicate=1, bamReads = "/Users/swethayadavalli/Desktop/SRR21677765/SRR21677765_sorted.bam")
Control_vs_Rescue<-dba.peakset(Control_vs_Rescue,
                          peaks="/Users/swethayadavalli/Desktop/SRR21677764_new_peaks_peaks.xls",
                          peak.caller="macs", sampID="Control-2",condition = "control", replicate=2, bamReads = "/Users/swethayadavalli/Desktop/SRR21677764/SRR21677764_sorted.bam")

Control_vs_Rescue<-dba.peakset(Control_vs_Rescue,
                          peaks="/Users/swethayadavalli/Desktop/SRR21677760_new_peaks_peaks.xls",
                          peak.caller="macs", sampID="Rescue-1",condition = "rescue", replicate=1, bamReads = "/Users/swethayadavalli/Desktop/SRR21677760/SRR21677760_sorted.bam")
Control_vs_Rescue<-dba.peakset(Control_vs_Rescue,
                          peaks="/Users/swethayadavalli/Desktop/SRR21677762_new_peaks_peaks.xls",
                          peak.caller="macs", sampID="Rescue-2",condition = "rescue", replicate=2, bamReads = "/Users/swethayadavalli/Desktop/SRR21677761/SRR21677761_sorted.bam")

Control_vs_Rescue_counts<-dba.count(Control_vs_Rescue, bParallel = TRUE, score=DBA_SCORE_READS)

png('/Users/swethayadavalli/Desktop/volcano.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(Control_vs_Rescue,Control_vs_Rescue$masks$control, main = "Open chromatic region overlaps in Control replicates")
dev.off()

png('/Users/swethayadavalli/Desktop/rescue.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(Control_vs_Rescue,Control_vs_Rescue$masks$rescue, main = "Open chromatic region overlaps in the Rescue replicates")
dev.off()


png('/Users/swethayadavalli/Desktop/heatmap.png', units="in", width=10, height=12, res=800)
plot(Control_vs_Rescue_counts, main="Correlation plot of studied samples")
dev.off()


Control_vs_Rescue_occupancy<-dba.peakset(Control_vs_Rescue, consensus = DBA_CONDITION,minOverlap = 0.33)
png('/Users/swethayadavalli/Desktop/occupancy_consensus_peaks_Control_vs_Rescue.png', units="in", width=7, heigh=7, res=800)
dba.plotVenn(Control_vs_Rescue_occupancy,Control_vs_Rescue_occupancy$masks$Consensus, main = "Overlap of consensus peaks of Control and Rescue")
dev.off()


#### Affinity analysis - differential accessibility analysis
counts<-dba.count(Control_vs_Rescue, bParallel = TRUE, score=DBA_SCORE_READS)
Control_vs_Rescue_counts<-dba.contrast(counts, categories=DBA_CONDITION,minMembers = 2)

Control_vs_Rescue_analyzed<-dba.analyze(Control_vs_Rescue_counts)
Control_vs_Rescue_DE_peaks <- dba.report(Control_vs_Rescue_analyzed)
Control_vs_Rescue_DE_peaks_DF<-as.data.frame(Control_vs_Rescue_DE_peaks)

write.table(Control_vs_Rescue_DE_peaks_DF[abs(Control_vs_Rescue_DE_peaks_DF$Fold)>1 & Control_vs_Rescue_DE_peaks_DF$FDR<0.05,], "/Users/swethayadavalli/Desktop/diff_regions_Rescue_vs_Control.txt", sep = "\t", quote = F)


