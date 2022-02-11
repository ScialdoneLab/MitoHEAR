#' get_raw_counts_allele
#'
#' It is one the two main function of the \strong{MitoHEAR} package (together
#' with \emph{get_heteroplasmy}). The function allows to obtain a matrix of
#' counts (n_row = number of sample, n_col= 4*number of bases) of the four
#' alleles in each base, for every sample. It takes as input a vector of sorted
#' bam files (one bam file for each sample) and a fasta file for the genomic
#' region of interest. It is based on the \emph{pileup} function of the package
#' Rsamtools.
#'
#'
#' @param bam_input Character vector of sorted bam files (full path). Each
#' sample is defined by one bam file. For each bam file it is needed also the
#' index bam file (.bai) at the same path.
#' @param path_fasta Character string with full path to the fasta file of the
#' genomic region of interest.
#' @param cell_names Character vector of sample names.
#' @param cores_number Number of cores to use.
#' @return A list with three elements:
#'
#' \item{matrix_allele_counts }{Matrix of counts (n_row = number of sample,
#' n_col= 4*number of bases) of the four alleles in each base, for every
#' sample. The row names is equal to cell_names. }
#' \item{name_position_allele}{Character vector with length equal to n_col of
#' matrix_allele_counts. Each element specifies the coordinate of genomic
#' position for a base and the allele.} \item{name_position}{Character vector
#' with length equal to n_col of matrix_allele_counts. Each element specifies
#' the coordinate of genomic position for a base.}
#' @author Gabriele Lubatti \email{gabriele.lubatti@@helmholtz-muenchen.de}
#' @seealso \url{https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/pileup}
#' @export get_raw_counts_allele
get_raw_counts_allele <- function(bam_input, path_fasta, cell_names, cores_number = 1) {


  fastaFile <- Biostrings::readDNAStringSet(path_fasta)
  seq_name <- names(fastaFile)
  sequence <- paste(fastaFile)
  df_bg <- data.frame(seq_name, sequence)
  name_sequences <- rep(0, length(df_bg$seq_name))
  for (i in 1:length(df_bg$seq_name)) {
    name_sequences[i] <- strsplit(df_bg$seq_name[i], split = " ")[[1]][1]
    message(paste0(name_sequences[i], " is in the provided fasta file"))
  }
  df_bg$seq_name <- name_sequences
  length_sequence <- rep(0, length(df_bg$seq_name))
  for (q in 1:length(df_bg$seq_name)) {
    length_sequence[q] <- length(unlist(strsplit(df_bg$sequence[q], split = NULL)))

    message(paste0(name_sequences[q], " has ", length_sequence[q], " bases"))
  }





  run_first_loop = function(t) {


    sbp <- Rsamtools::ScanBamParam(which = GenomicRanges::GRanges(name_sequences, IRanges::IRanges(rep(1, length(df_bg$seq_name)), length_sequence)))

    chr_samtools <- Rsamtools::pileup((bam_input[t]), scanBamParam = sbp, pileupParam = Rsamtools::PileupParam(min_nucleotide_depth = 0, min_base_quality =  33, distinguish_strands = FALSE, include_insertions = FALSE, distinguish_nucleotides = TRUE, ignore_query_Ns = FALSE, include_deletions = FALSE))


    counts_samtools <- reshape2::dcast(chr_samtools, seqnames + pos ~ nucleotide, value.var = "count", fun.aggregate = sum)

    chr_names <- df_bg$seq_name
    list_chr_all <- as.list(rep(0, length(chr_names[chr_names%in%as.vector(counts_samtools$seqnames)]) ))

    for ( k in 1:length(chr_names[chr_names%in%as.vector(counts_samtools$seqnames)])) {
      chr_selected <- chr_names[chr_names%in%as.vector(counts_samtools$seqnames)][k]
      counts_samtools_chr <- counts_samtools[grep(chr_selected, as.vector(counts_samtools$seqnames)), ]
      chr_seq <- unlist(strsplit(df_bg$sequence[df_bg$seq_name==chr_selected], split = NULL))
      base_spikes <- chr_seq
      allele_1 <- rep(0, length(chr_seq))
      allele_2 <- rep(0, length(chr_seq))
      allele_3 <- rep(0, length(chr_seq))
      allele_4 <- rep(0, length(chr_seq))
      bases_chr <- seq(1:length(chr_seq))
        for (i in counts_samtools_chr$pos) {
          allele_1[i] <- counts_samtools_chr$A[counts_samtools_chr$pos==i]
          allele_2[i] <- counts_samtools_chr$C[counts_samtools_chr$pos==i]
          allele_3[i] <- counts_samtools_chr$G[counts_samtools_chr$pos==i]
          allele_4[i] <- counts_samtools_chr$T[counts_samtools_chr$pos==i]
        }

      final_counts_samtools <- data.frame(rep(chr_selected, length(chr_seq)), bases_chr, allele_1, allele_2, allele_3, allele_4, base_spikes)
      colnames(final_counts_samtools) <- c("seqnames", "pos", "A", "C", "G", "T", "Reference")


      name_A <- paste(final_counts_samtools$pos, "A", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")
      name_B <- paste(final_counts_samtools$pos, "C", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")
      name_C <- paste(final_counts_samtools$pos, "G", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")
      name_D <- paste(final_counts_samtools$pos, "T", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")

      chr_seq <- unlist(strsplit(df_bg$sequence[df_bg$seq_name==chr_selected], split=NULL))
      select_allele <- apply(final_counts_samtools, 1, function(x) {

        x <- x[3:6]
        x <- as.numeric(x)
        return(x)
      })

      allele_matrix <- matrix(0, ncol = 4*length(row.names(final_counts_samtools)), nrow = 1)
      j <- 1
      for (i in 1:dim(select_allele)[2]) {
        allele_matrix[, j:(j + 3)] = select_allele[, i]

        j <- j + 4
      }

      j <- 1
      name_allele <- rep(0, length(colnames(allele_matrix)))
      for (i in 1:length(row.names(final_counts_samtools))) {
        name_allele[j:(j +  3)] = c(name_A[i], name_B[i], name_C[i], name_D[i])
        j <- j +  4
      }

      colnames(allele_matrix) <- name_allele
      row.names(allele_matrix) <- cell_names[t]
      list_chr_all[[k]] <- allele_matrix
      names(list_chr_all[[k]]) <- chr_selected



    }

    message(paste0("Sample ", cell_names[t]))

    return(list_chr_all)

  }





  list_all_cell <- mclapply(seq(1:length(bam_input)), run_first_loop ,  mc.cores = cores_number)



  list_chr_names <- as.list(rep(0, length(bam_input)))

  for (t in 1:length(bam_input)) {
    name_chr <- rep(0, length(list_all_cell[[t]]))
    for (k in 1:length(name_chr)) {


      name_chr[k] <- names(list_all_cell[[t]][[k]][1])}
    list_chr_names[[t]] <- name_chr
  }


  all_chr <- unlist(list_chr_names)
  all_chr <- unique(all_chr)





  run_second_loop = function(t) {

    absent_chr <- all_chr[!all_chr%in%list_chr_names[[t]]]

    list_chr_absent <- as.list(rep(0, length(absent_chr)))
    if (length(absent_chr)>0) {


      for ( k in 1:length(absent_chr)) {


        chr_selected <- absent_chr[k]

        warning(paste0("Chr ", absent_chr[k], "is not present in sample: ", cell_names[t]))
        chr_seq <- unlist(strsplit(df_bg$sequence[df_bg$seq_name==chr_selected], split = NULL))

        allele_1 <- rep(0, length(chr_seq))
        allele_2 <- rep(0, length(chr_seq))
        allele_3 <- rep(0, length(chr_seq))
        allele_4 <- rep(0, length(chr_seq))
        bases_chr <- seq(1:length(chr_seq))

        final_counts_samtools <- data.frame(rep(chr_selected, length(chr_seq)), bases_chr, allele_1, allele_2, allele_3, allele_4, chr_seq)
        colnames(final_counts_samtools) <- c("seqnames", "pos", "A", "C", "G", "T", "Reference")




        name_A <- paste(final_counts_samtools$pos, "A", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")
        name_B <- paste(final_counts_samtools$pos, "C", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")
        name_C <- paste(final_counts_samtools$pos, "G", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")
        name_D <- paste(final_counts_samtools$pos, "T", final_counts_samtools$Reference, final_counts_samtools$seqnames, sep = "_")


        select_allele <- apply(final_counts_samtools, 1, function(x) {

          x <- x[3:6]
          x <- as.numeric(x)
          return(x)
        })

        allele_matrix <- matrix(0, ncol = 4*length(row.names(final_counts_samtools)), nrow = 1)
        j <- 1
        for (i in 1:dim(select_allele)[2]) {
          allele_matrix[, j:(j + 3)] <- select_allele[, i]

          j <- j + 4
        }

        j <- 1
        name_allele <- rep(0, length(colnames(allele_matrix)))
        for (i in 1:length(row.names(final_counts_samtools))) {
          name_allele[j:(j + 3)] <- c(name_A[i], name_B[i], name_C[i], name_D[i])
          j <- j + 4
        }

        colnames(allele_matrix) <- name_allele
        row.names(allele_matrix) <- cell_names[t]
        list_chr_absent[[k]] <- allele_matrix


      }
    }

    else {}



    return(list_chr_absent)

  }



  list_all_cell_absent <- mclapply(seq(1:length(bam_input)), run_second_loop ,  mc.cores = cores_number)





  run_third_loop = function(t) {



    list_complete <- c(list_all_cell[[t]], list_all_cell_absent[[t]])
    return(list_complete)
  }
  list_all_complete <- mclapply(seq(1:length(bam_input)), run_third_loop ,  mc.cores = cores_number)


  chr_length <- rep(0, length(all_chr))
  for (t in 1:length(all_chr)) {
    chr_length[t] <- length(unlist(strsplit(df_bg$sequence[df_bg$seq_name==all_chr[t]], split = NULL)))
  }



  matrix_allele_counts <- matrix(0, ncol = sum(chr_length)*4, nrow = length(bam_input))

  for (i in 1:length(bam_input)) {

    matrix_allele_counts[i, ] <- rlist::list.cbind(list_all_complete[[i]])
  }

  colnames(matrix_allele_counts) <- colnames(rlist::list.cbind(list_all_complete[[1]]))
  row.names(matrix_allele_counts) <- cell_names


  name_position_allele <- colnames(matrix_allele_counts)
  name_position <- rep(0, length(colnames(matrix_allele_counts)))
  for (i in 1:length(colnames(matrix_allele_counts))) {
    part_1 <- strsplit(colnames(matrix_allele_counts)[i], split = "_")[[1]][1]
    part_2 <- strsplit(colnames(matrix_allele_counts)[i], split = "_")[[1]][4]
    name_position[i] <- paste(part_1, part_2, sep = "_")

  }


  name_position_unique <- unique(name_position)



  return(list(matrix_allele_counts, name_position_allele, name_position))

}


