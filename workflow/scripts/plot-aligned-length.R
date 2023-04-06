#!/usr/bin/env Rscript
#
# Visualize length of mapped reads (total) and length of aligned portion of reads (excludes softclipping, etc)
#
# Four arguments: length of mapped reads, length of aligned portion of reads, read length cap for the plot, output pdf file
#

suppressPackageStartupMessages(library("rio"))

#args <- commandArgs(trailingOnly = TRUE)

# args <- NULL
# args <- c("/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.mapped-len.full.tsv.gz", # Whole mapped read length
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.align-len.full.tsv.gz", # Only the aligned portion
#           5000, # Maximum read length for plotting
#           "/home/jan/playground/TERA-Seq_snakemake/data/samples/hsa.dRNASeq.HeLa.polyA.CIP.decap.REL5.long.1/log/reads.1.sanitize.toGenome.len.mapped-vs-aligned.pdf")

mapped_len <- snakemake@input[["mapped_len"]] # args[1]
aligned_portion_len <- snakemake@input[["aligned_len"]] #  args[2]
read_cap <- snakemake@params[["len_cap"]] #  as.numeric(args[3])
ofile <- snakemake@output[[1]] #  args[4]

read_cap <- as.numeric(read_cap)

aligned_portion_len <- rio::import(aligned_portion_len, format = "tsv")
mapped_len <- rio::import(mapped_len, format = "tsv")

if(nrow(aligned_portion_len) == 0 | nrow(mapped_len) == 0){
  pdf(ofile)
    par(mar = c(0, 0, 0, 0))
    plot(c(0, 1), c(0, 1), ann = F, bty = "n", type = "n", xaxt = "n", yaxt = "n")
    text(
      x = 0.5, y = 0.5, paste("Either\n", snakemake@input[["mapped_len"]], "\nor\n", snakemake@input[["aligned_len"]], "\ndidn't contain any lines. Please check the input files."),
      cex = 1, col = "black"
    )
  dev.off()
} else {
  aligned_portion_len <- aligned_portion_len[aligned_portion_len$length != 1, ] # Reads with length 1 are usually secondary alignments which have "*" instead of read sequence
  
  # Get median length
  med_tot <- median(mapped_len$length)
  med_mapped <- median(aligned_portion_len$length)
  
  aligned_portion_len <- as.matrix(table(aligned_portion_len$length))
  mapped_len <- as.matrix(table(mapped_len$length))
  
#  xmax <- max(c(as.numeric(rownames(mapped_len)), as.numeric(rownames(aligned_portion_len))))
  xmax <- read_cap
  ymax<- max(c(mapped_len, aligned_portion_len))

  pdf(ofile)
    plot(
      x = as.numeric(rownames(mapped_len)), y = mapped_len, cex = 0.2, col = "red", type = "l",
      xlim = c(0, xmax), ylim = c(0, ymax * 1.05),
      xlab = paste0("Length (capped at ", xmax, " nt)"), ylab = "Number of reads",
      main = paste0("Mapped read vs Aligned portion length\n", "\nMedian length - Mapped: ", med_tot, "; Aligned portion: ", med_mapped)
    )
    lines(x = as.numeric(rownames(aligned_portion_len)), y = aligned_portion_len, cex = 0.2, col = "blue")
    legend("topright",
      legend = c("Mapped", "Aligned portion"),
      col = c("red", "blue"), lty = 1, cex = 1.2
    )
  dev.off()
}
