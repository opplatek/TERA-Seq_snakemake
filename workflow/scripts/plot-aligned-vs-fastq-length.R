#!/usr/bin/env Rscript
#
# Compare length of all reads vs mapped reads (binned by 50 nt) and the first bin to reach 90% of mapped reads
#
# Note: This will not work very well for datasets with low number of reads
# Note: "First to get 90% aligned length" is calculated from binned lengths (by 50 nt)
#
# Four arguments: input total length of trimmed (or raw) reads, only aligned read length, read length cap, output pdf
#

library("rio")
suppressPackageStartupMessages(library("dplyr"))

args <- commandArgs(trailingOnly = TRUE)

# args <- NULL
# args <- c("/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.adapt_trim.read-len.tsv.gz", # Trimme fastq read length (incl. unmapped)
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz", # Whole mapped read length
#           5000,  # Maximum read length for plotting
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.len.trim-vs-mapped.pdf")

total_len <- args[1]
aligned_len <- args[2]
read_cap <- as.numeric(args[3])
ofile <- args[4]

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
  
  colnames(total_len)[colnames(total_len) == "length"] <- "total_length"
  colnames(aligned_len)[colnames(aligned_len) == "length"] <- "aln_length"
  
  len <- merge(total_len, aligned_len, all.x = T, by = "read_id")
  len[is.na(len)] <- 0
  
  # Get median length
  med_tot <- median(len$total_length, na.rm = T) # Expand the frequency tab
  med_aln <- median(len$aln_length, na.rm = T) # Expand the frequency tab
  
  # Make bins by 50 nucleotides
  bin_len <- 50
  
  len <- len %>%
    filter(total_length <= !!read_cap)
  
  len$ratio <- len$aln_length / (len$total_length / 100)
  
  len_binned <- len %>%
    group_by(group = cut(total_length, breaks = seq(0, max(total_length), !!bin_len), labels = seq(!!bin_len, max(total_length), !!bin_len))) %>%
    summarize("ratio"=median(ratio, na.rm = T)) %>%
    rename("length" = group) %>%
    as.data.frame()
  
  first_ninety <- as.numeric(as.character(len_binned[len_binned$ratio >= 90, "length"][1]))
  
  len <- len %>%
    group_by(total_length) %>%
    summarize("ratio" = median(ratio))
  
  pdf(ofile)
    plot(
      x = len$total_length, y = len$ratio,
      axes = F,
      cex = 0.2, col = "red",
      xlab = "Read length", ylab = "Perc. of aligned",
      main = paste0("Mapped vs All read length ratios\n", "Median length - All: ", med_tot, "; Mapped total: ", med_aln, "\nFirst len. to have >=90% mapped: ", paste0(c(first_ninety-bin_len, first_ninety), collapse = "-"), " nt")
    )
  #  axis(1, at = seq(0, xmax, by = 250), labels = seq(0, xmax, by = 250))
    axis(1, at = seq(0, read_cap, by=200), labels = seq(0, read_cap, by=200))
    axis(2, at = seq(0, 100, by = 10), labels = seq(0, 100, by = 10))
    abline(v = c(first_ninety-bin_len, first_ninety), col = "blue", lty = 3) # First length to reach 90% aln. rate
    abline(h = mean(len$ratio), col = "black", lty = 2)
    abline(h = median(len$ratio), col = "black")
    legend("center",
      legend = c("Mapped/All ratio", "Ratio mean", "Ratio median", ">=90% mapped"),
      col = c("red", "black", "black", "blue"), pch = c(20, NA, NA, NA), lty = c(NA, 2, 1, 3), lwd = c(NA, 1, 1, 1)
    )
  dev.off()
}