fix_peak_range <- function(peak_info, target_width){
  l_width <- peak_info$summit - peak_info$start
  r_width <- peak_info$end - peak_info$summit
  new_start <- round(peak_info$summit - target_width * l_width / (l_width + r_width))
  new_end <- round(peak_info$summit + target_width * r_width / (l_width + r_width))
  peak_info$start <- new_start
  peak_info$end <- new_end
  peak_info
}


get_peaks_altius <- function(altius_peak_mat_path, altius_sample_info_path, altius_peak_info_path, n_systems, fixed_peak_width = NULL){
  require("Matrix")
  require("MatrixGenerics")
  require("GenomicRanges")
  
  #files from here: https://zenodo.org/records/3838751
  altius_peak_mat <- Matrix::readMM(altius_peak_mat_path)
  #for some reason, line 734 is malformed? It shouldn't be there anyway
  altius_sample_info <- read.delim(altius_sample_info_path, sep = "\t", header = T)[1:733,]
  altius_peak_info <- read.delim(altius_peak_info_path, sep = "\t", header = T)
  
  mat_list <- list()
  for(organ.system in sort(unique(altius_sample_info$System))){
    mat_list[[organ.system]] <- altius_peak_mat[,altius_sample_info$System == organ.system]
  }
  system_mat <- do.call(cbind, lapply(mat_list, rowSums))
  use_idx <- which(rowSums(system_mat != 0) == n_systems)
  
  altius_peak_info <- altius_peak_info[use_idx,]
  if (!is.null(fixed_peak_width) & is.integer(fixed_peak_width)){
    altius_peak_info <- fix_peak_range(altius_peak_info, fixed_peak_width)
  }
  GenomicRanges::makeGRangesFromDataFrame(altius_peak_info)
}
