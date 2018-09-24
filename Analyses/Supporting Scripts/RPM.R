RPM <- function(counts){
  #counts is matrix with samples in columns and genes in rows
  total_reads <- colSums(counts)
  rpm <- t(apply(counts, 2, function(x) x*1e6/sum(x)))
  detected_trans <- rowSums(rpm>1)
  stuff <- list("total_reads"= total_reads, "RPM" = rpm, "log2RPM" = log2(1+rpm), "detected_trans"= detected_trans)
  return(stuff)
}