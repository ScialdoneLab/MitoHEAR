# Example with two samples and two bases whose reference allele is A and C.
# The two samples have 100 reads in the reference allele and 0 in all the others.

sample1_A <- c(1, 0, 0, 0)
sample1_C <- c(1, 0, 0, 0)
sample1 <- c(sample1_A, sample1_C)
sample2_A <- c(1, 0, 0, 0)
sample2_C <- c(1, 0, 0, 0)
sample2 <- c(sample2_A, sample2_C)


names_A <- rep("1_A", length(sample1_A))
names_C <- rep("2_C", length(sample1_C))
allele <- c("A", "C", "T", "G")
names_A_allele <- paste(names_A, allele, sep = " ")
names_C_allele <- paste(names_C, allele, sep = " ")

test_allele <- matrix(c(sample1, sample2), byrow = T, ncol = 8, nrow = 2)
colnames(test_allele) <- c(names_A_allele, names_C_allele)
row.names(test_allele) <- c("sample1", "sample2")


name_position_allele_test <- c(names_A_allele, names_C_allele)
name_position_test <- c(names_A, names_C)

cluster_test <- c("cluster_test", "cluster_test")

# Build heteroplasmy matrix for the 2 samples.
heteroplasmy_matrix_test <- matrix(0, ncol = length(unique(name_position_test)), nrow = length(row.names(test_allele)))
colnames(heteroplasmy_matrix_test) <- unique(name_position_test)
row.names(heteroplasmy_matrix_test) <- c("sample1", "sample2")

# Build allele matrix for the 2 samples.
allele_matrix_test <- matrix(c(sample1, sample2), byrow = T, ncol = length(name_position_test), nrow = length(row.names(test_allele)))
colnames(allele_matrix_test) <- c(names_A_allele, names_C_allele)
row.names(allele_matrix_test) <- c("sample1", "sample2")



test_that("clustering_angular_distance gives error when samples have the same allele frequencies matrix",  {
  expect_error(clustering_angular_distance(heteroplasmy_matrix_test, allele_matrix_test, cluster_test, length(row.names(heteroplasmy_matrix_test)), deepSplit_param = 0,
                                           minClusterSize_param = 1, 0.2, min_value = 0.001, index = NULL, relevant_bases = NULL), "All the bases has variance 0 across samples")

})


