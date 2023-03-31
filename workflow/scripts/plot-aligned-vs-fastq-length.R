#!/usr/bin/env Rscript
#
# Compare length of all reads vs mapped reads (binned by 50 nt) and the first bin to reach 90% of mapped reads
#
# Note: This will not work very well for datasets with low number of reads
#

library("rio")
suppressPackageStartupMessages(library("dplyr"))

args <- commandArgs(trailingOnly = TRUE)

# args <- NULL
# args <- c("/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.adapt_trim.read-len.txt", # Trimme fastq read length (incl. unmapped)
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.align-len.full.txt", # Whole mapped read length
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.len.trim-vs-mapped.pdf")

total_len <- args[1]
aligned_len <- args[2]
ofile <- args[3]

total_len <- rio::import(total_len, format = "tsv")
aligned_len <- rio::import(aligned_len, format = "tsv")

if(nrow(total_len) == 0 | nrow(aligned_len) == 0){
  pdf(ofile)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    text(
      x = 0.5, y = 0.5, paste("Either\n", args[1], "\nor\n", args[2], "\ndidn't contain any lines. Please check the input files."),
      cex = 1, col = "black"
    )
  dev.off()
} else {
  aligned_len <- aligned_len[aligned_len$length != 1, ] # Reads with length 1 are usually secondary alignments which have "*" instead of read sequence
  
  aligned_len <- as.data.frame(table(aligned_len$length)) # Make frequency table
  colnames(aligned_len) <- c("length", "count")
  
  colnames(total_len)[colnames(total_len) == "count"] <- "total_ct"
  colnames(aligned_len)[colnames(aligned_len) == "count"] <- "aln_ct"
  
  len <- merge(total_len, aligned_len, all.x = T, by = "length")
  len[is.na(len)] <- 0
  
  # Get median length
  med_tot <- median(rep(len$length, len$total_ct)) # Expand the frequency tab
  med_aln <- median(rep(len$length, len$aln_ct)) # Expand the frequency tab
  
  # Make bins by 50 nucleotides
  bin_len <- 50
  
  len <- len %>%
    group_by(group = cut(length, breaks = seq(0, max(length), !!bin_len))) %>%
    summarise(total_ct = sum(total_ct), aln_ct = sum(aln_ct)) %>%
    rename("length" = group) %>%
    as.data.frame()
  
  len$ratio <- len$aln_ct / (len$total_ct / 100)
  
  pdf(ofile)
    plot(
      x = len$length, y = len$ratio,
      axes = F,
      cex = 0.2, col = "red",
      xlab = "Read length", ylab = "Perc. of aligned",
      main = paste0("Mapped vs All read length ratios\n", "Median length - All: ", med_tot, "; Mapped total: ", med_aln, "\nFirst len. to have >=90% mapped: ", len[len$ratio >= 90, "length"][1])
    )
  #  axis(1, at = seq(0, xmax, by = 250), labels = seq(0, xmax, by = 250))
    axis(1, at = len$length, labels = len$length)
    axis(2, at = seq(0, 100, by = 10), labels = seq(0, 100, by = 10))
    abline(v = len[len$ratio > 90, "length"][1], col = "blue", lty = 3) # First length to reach 90% aln. rate
    abline(h = mean(len$ratio), col = "black", lty = 2)
    abline(h = median(len$ratio), col = "black")
    legend("topright",
      legend = c("Mapped/All ratio", "Ratio mean", "Ratio median", ">=90% mapped"),
      col = c("black", "black", "black", "blue"), pch = c(NA, NA, NA, NA), lty = c(1, 2, 1, 3), lwd = c(5, 1, 1, 1)
    )
  dev.off()
}