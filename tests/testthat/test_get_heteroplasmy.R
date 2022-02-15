# Example with two samples and two bases whose reference allele is A and C.
# The two samples have 100 reads in the reference allele and 0 in all the others.

sample1_A=c(100,0,0,0)
names_A=rep("1_A",length(sample1_A))
sample1_C=c(100,0,0,0)
names_C=rep("2_C",length(sample1_C))
allele=c("A","C","T","G")
names_A_allele=paste(names_A,allele,sep=" ")
names_C_allele=paste(names_C,allele,sep=" ")
sample1=c(sample1_A,sample1_C)
sample2_A=c(100,0,0,0)
sample2_C=c(100,0,0,0)
sample2=c(sample2_A,sample2_C)
test_allele=matrix(c(sample1,sample2),byrow = T,ncol=8,nrow=2)
colnames(test_allele)=c(names_A_allele,names_C_allele)
row.names(test_allele)=c("sample1","sample2")

name_position_allele_test=c(names_A_allele,names_C_allele)
name_position_test=c(names_A,names_C)

# Build the list of output
test_output_1=matrix(100,ncol=length(unique(name_position_test)),nrow=length(row.names(test_allele)))
colnames(test_output_1)=unique(name_position_test)
row.names(test_output_1)=c("sample1","sample2")

test_output_2=matrix(100,ncol=length(unique(name_position_test)),nrow=length(row.names(test_allele)))
colnames(test_output_2)=unique(name_position_test)
row.names(test_output_2)=c("sample1","sample2")

test_output_3=matrix(0,ncol=length(unique(name_position_test)),nrow=length(row.names(test_allele)))
colnames(test_output_3)=unique(name_position_test)
row.names(test_output_3)=c("sample1","sample2")

sample1_A=c(1,0,0,0)
sample1_C=c(1,0,0,0)
sample1=c(sample1_A,sample1_C)
sample2_A=c(1,0,0,0)
sample2_C=c(1,0,0,0)
sample2=c(sample2_A,sample2_C)
test_output_4=matrix(c(sample1,sample2),byrow = T,ncol=length(name_position_test),nrow=length(row.names(test_allele)))
colnames(test_output_4)=c(names_A_allele,names_C_allele)
row.names(test_output_4)=c("sample1","sample2")

test_output_5=NULL
test_output=list(test_output_1,test_output_2,test_output_3,test_output_4,test_output_5)


test_that("get_heteroplasmy works with a small dataset of two samples and two bases", {
  expect_equal(get_heteroplasmy(test_allele,name_position_allele_test,name_position_test,number_reads=50,number_positions=1,filtering = 1),test_output)
})

test_that("get_heteroplasmy gives error when there are not at least 2 samples that cover at least 2 bases with more than 100 reads", {
  expect_error(get_heteroplasmy(test_allele,name_position_allele_test,name_position_test,number_reads=100,number_positions=1,filtering = 1))
})

test_that("get_heteroplasmy gives error when there are not at least 2 samples and at least 2 bases as input", {
  expect_error(get_heteroplasmy(test_allele[1,],names_A_allele,names_A,number_reads=50,number_positions=1,filtering = 1),paste0('There are not at least 2 samples and at least 2 bases as input'))
})


