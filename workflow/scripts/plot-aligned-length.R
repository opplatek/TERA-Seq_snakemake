#!/usr/bin/env Rscript
#
# Visualize length of mapped reads (total) and length of mapped portion of reads (excludes softclipping, etc)
#
# Three arguments: length of mapped reads, length of mapped portion of reads, output pdf file
#

library("rio")

args <- commandArgs(trailingOnly = TRUE)

# args <- NULL
# args <- c("/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.align-len.full.txt", # Whole mapped read length
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.mapped-len.full.txt", # Only the aligned portion
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.len.mapped.pdf")

aligned_len <- args[1]
mapped_portion_len <- args[2]
ofile <- args[3]

mapped_portion_len <- rio::import(mapped_portion_len, format = "tsv")
aligned_len <- rio::import(aligned_len, format = "tsv")

if(nrow(mapped_portion_len) == 0 | nrow(aligned_len) == 0){
  pdf(ofile)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    text(
      x = 0.5, y = 0.5, paste("Either\n", args[1], "\nor\n", args[2], "\ndidn't contain any lines. Please check the input files."),
      cex = 1, col = "black"
    )
  dev.off()
} else {
  mapped_portion_len <- mapped_portion_len[mapped_portion_len$length != 1, ] # Reads with length 1 are usually secondary alignments which have "*" instead of read sequence
  
  # Get median length
  med_tot <- median(aligned_len$length)
  med_mapped <- median(mapped_portion_len$length)
  
  mapped_portion_len <- as.matrix(table(mapped_portion_len$length))
  aligned_len <- as.matrix(table(aligned_len$length))
  
  xmax <- max(as.numeric(rownames(aligned_len)), as.numeric(rownames(mapped_portion_len)))
  
  pdf(ofile)
    plot(
      x = as.numeric(rownames(aligned_len)), y = aligned_len, cex = 0.2, col = "red", type = "l",
      xlim = c(0, xmax),
      xlab = paste0("Length (capped at ", xmax, " nt)"), ylab = "Number of reads",
      main = paste0("Aligned portion vs Mapped read length\n", "\nMedian length - Aligned portion: ", med_mapped,  "; Mapped: ", med_tot)
    )
    lines(x = as.numeric(rownames(mapped_portion_len)), y = mapped_portion_len, cex = 0.2, col = "blue")
    legend("topright",
      legend = c("Aligned portion", "Mapped"),
      col = c("blue", "red"), lty = 1, cex = 1.2
    )
  dev.off()
}
