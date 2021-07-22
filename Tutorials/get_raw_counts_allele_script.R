
gc()
rm(list=ls())
library(optparse)
option_list = list(
  make_option(c("-b", "--bam_input"), type="character", default=NULL,
              help="character vector of sorted bam files", metavar="character"),
  make_option(c("-f", "--path_fasta"), type="character", default=NULL,
              help="fasta file of the genomic region of interest", metavar="integer"),
  make_option(c("-c", "--cell_names"), type="character", default=NULL,
              help="character vector with the names of the samples that will be used in the down-stream analysis", metavar="character"),
  make_option(c("-o", "--output"), type="character", default="out.Rda",
              help="output file name [default= %default]", metavar="character"),
  make_option(c("-s", "--cores_number"), type="integer", default=1,
              help="Number of cores to use[default= %default]", metavar="integer")
)

opt_parser = OptionParser(option_list=option_list,add_help_option=FALSE)
opt = parse_args(opt_parser)

if (is.null(opt$'cell_names')|is.null(opt$'bam_input')| is.null(opt$'path_fasta')){
  print_help(opt_parser)
  stop("Three arguments must be supplied (bam_input,path_fasta,cell_names).n", call.=FALSE)
}

load(paste0(opt$'cell_names',".Rda"))
load(paste0(opt$'bam_input',".Rda"))






library(Rsamtools)







get_raw_counts_allele=function(bam_input,path_fasta,cell_names,cores_number=1){


  fastaFile <- Biostrings::readDNAStringSet(path_fasta)
  seq_name = names(fastaFile)
  sequence = paste(fastaFile)
  df_bg <- data.frame(seq_name, sequence)
  name_sequences=rep(0,length(df_bg$seq_name))
  for (i in 1:length(df_bg$seq_name)){
    name_sequences[i]=strsplit(df_bg$seq_name[i],split=" ")[[1]][1]
    print(name_sequences[i])
  }
  df_bg$seq_name=name_sequences
  length_sequence=rep(0,length(df_bg$seq_name))
  for (q in 1:length(df_bg$seq_name)){
    length_sequence[q]=length(unlist(strsplit(df_bg$sequence[q],split=NULL)))
    print(length_sequence[q])
  }





  run_first_loop=function(t){


    sbp <- Rsamtools::ScanBamParam(which=GRanges(name_sequences,IRanges(rep(1,length(df_bg$seq_name)),length_sequence)))

    example_spikes=Rsamtools::pileup((bam_input[t]),scanBamParam=sbp,pileupParam=PileupParam(min_nucleotide_depth=0,min_base_quality= 33, distinguish_strands=FALSE,include_insertions=FALSE,distinguish_nucleotides=TRUE,ignore_query_Ns=FALSE,include_deletions=FALSE))


    matrix_example_spikes_2=data.table::dcast(example_spikes, seqnames+pos ~ nucleotide, value.var="count", fun.aggregate=sum)

    ercc_name_full=df_bg$seq_name
    list_all_allele=as.list(rep(0,length(ercc_name_full[ercc_name_full%in%as.vector(matrix_example_spikes_2$seqnames)]) ))

    for ( k in 1:length(ercc_name_full[ercc_name_full%in%as.vector(matrix_example_spikes_2$seqnames)])){


      ercc_name=ercc_name_full[ercc_name_full%in%as.vector(matrix_example_spikes_2$seqnames)][k]
      new_matrix=matrix_example_spikes_2[grep(ercc_name,as.vector(matrix_example_spikes_2$seqnames)),]
      spikes_1=unlist(strsplit(df_bg$sequence[df_bg$seq_name==ercc_name],split=NULL))
      base_spikes=spikes_1
      allele_1=rep(0,length(spikes_1))
      allele_2=rep(0,length(spikes_1))
      allele_3=rep(0,length(spikes_1))
      allele_4=rep(0,length(spikes_1))
      pos_new=seq(1:length(spikes_1))
      for (i in new_matrix$pos){
        allele_1[i]=new_matrix$A[new_matrix$pos==i]
        allele_2[i]=new_matrix$C[new_matrix$pos==i]
        allele_3[i]=new_matrix$G[new_matrix$pos==i]
        allele_4[i]=new_matrix$T[new_matrix$pos==i]

      }

      final_matrix=data.frame(rep(ercc_name,length(spikes_1)),pos_new,allele_1,allele_2,allele_3,allele_4,base_spikes)
      colnames(final_matrix)=c("seqnames","pos","A","C","G","T","Reference")

      nuova_matrice=matrix(0,ncol = 4*length(row.names(new_matrix)),nrow = 1)


      name_A=paste(final_matrix$pos,"A",final_matrix$Reference,final_matrix$seqnames,sep="_")
      name_B=paste(final_matrix$pos,"C",final_matrix$Reference,final_matrix$seqnames,sep="_")
      name_C=paste(final_matrix$pos,"G",final_matrix$Reference,final_matrix$seqnames,sep="_")
      name_D=paste(final_matrix$pos,"T",final_matrix$Reference,final_matrix$seqnames,sep="_")

      spikes_1=unlist(strsplit(df_bg$sequence[df_bg$seq_name==ercc_name],split=NULL))
      so=apply(final_matrix,1,function(x){

        x=x[3:6]
        x=as.numeric(x)
        return(x)
      })

      allele_matrix=matrix(0,ncol = 4*length(row.names(final_matrix)),nrow = 1)
      j=1
      for (i in 1:dim(so)[2]){
        allele_matrix[,j:(j+3)]=so[,i]

        j=j+4}

      j=1
      name_full=rep(0,length(colnames(allele_matrix)))
      for (i in 1:length(row.names(final_matrix))){
        name_full[j:(j+3)]=c(name_A[i],name_B[i],name_C[i],name_D[i])
        j=j+4
      }

      colnames(allele_matrix)=name_full
      row.names(allele_matrix)=cell_names[t]
      list_all_allele[[k]]=allele_matrix
      names(list_all_allele[[k]])=ercc_name

      print(paste(k,ercc_name,sep=""))

    }
    print(paste(t,cell_names[t]))

    return(list_all_allele)

  }





  list_all_cell=parallel::mclapply(seq(1:length(bam_input)),run_first_loop , mc.cores = cores_number)



  list_all_cell_spikes=as.list(rep(0,length(bam_input)))

  for (t in 1:length(bam_input)){
    name_spikes_cell=rep(0,length(list_all_cell[[t]]))
    for (k in 1:length(name_spikes_cell)){


      name_spikes_cell[k]=names(list_all_cell[[t]][[k]][1])}
    list_all_cell_spikes[[t]]=name_spikes_cell
  }


  all_spikes=unlist(list_all_cell_spikes)
  all_spikes=unique(all_spikes)





  run_second_loop=function(t){


    absent_spikes=all_spikes[!all_spikes%in%list_all_cell_spikes[[t]]]

    list_all_allele=as.list(rep(0,length(absent_spikes)))
    if (length(absent_spikes)>0){
      print(absent_spikes)

      for ( k in 1:length(absent_spikes)){


        ercc_name=absent_spikes[k]
        print(ercc_name)
        spikes_1=unlist(strsplit(df_bg$sequence[df_bg$seq_name==ercc_name],split=NULL))
        base_spikes=spikes_1
        allele_1=rep(0,length(spikes_1))
        allele_2=rep(0,length(spikes_1))
        allele_3=rep(0,length(spikes_1))
        allele_4=rep(0,length(spikes_1))
        pos_new=seq(1:length(spikes_1))

        final_matrix=data.frame(rep(ercc_name,length(spikes_1)),pos_new,allele_1,allele_2,allele_3,allele_4,base_spikes)
        colnames(final_matrix)=c("seqnames","pos","A","C","G","T","Reference")

        nuova_matrice=matrix(0,ncol = 4*length(row.names(final_matrix)),nrow = 1)


        name_A=paste(final_matrix$pos,"A",final_matrix$Reference,final_matrix$seqnames,sep="_")
        name_B=paste(final_matrix$pos,"C",final_matrix$Reference,final_matrix$seqnames,sep="_")
        name_C=paste(final_matrix$pos,"G",final_matrix$Reference,final_matrix$seqnames,sep="_")
        name_D=paste(final_matrix$pos,"T",final_matrix$Reference,final_matrix$seqnames,sep="_")

        spikes_1=unlist(strsplit(df_bg$sequence[df_bg$seq_name==ercc_name],split=NULL))
        so=apply(final_matrix,1,function(x){

          x=x[3:6]
          x=as.numeric(x)
          return(x)
        })

        allele_matrix=matrix(0,ncol = 4*length(row.names(final_matrix)),nrow = 1)
        j=1
        for (i in 1:dim(so)[2]){
          allele_matrix[,j:(j+3)]=so[,i]

          j=j+4}

        j=1
        name_full=rep(0,length(colnames(allele_matrix)))
        for (i in 1:length(row.names(final_matrix))){
          name_full[j:(j+3)]=c(name_A[i],name_B[i],name_C[i],name_D[i])
          j=j+4
        }

        colnames(allele_matrix)=name_full
        row.names(allele_matrix)=cell_names[t]
        list_all_allele[[k]]=allele_matrix


      }}
    else{}

    print(paste(t,cell_names[t]))

    return(list_all_allele)

  }



  list_all_cell_absent=parallel::mclapply(seq(1:length(bam_input)),run_second_loop , mc.cores = cores_number)





  run_third_loop=function(t){



    list_all_complete_one=c(list_all_cell[[t]],list_all_cell_absent[[t]])
    return(list_all_complete_one)
    }
  list_all_complete=parallel::mclapply(seq(1:length(bam_input)),run_third_loop , mc.cores = cores_number)


  spikes_length=rep(0,length(all_spikes))
  for (t in 1:length(all_spikes)){
    spikes_length[t]=length(unlist(strsplit(df_bg$sequence[df_bg$seq_name==all_spikes[t]],split=NULL)))
  }



  matrix_allele_counts=matrix(0,ncol = sum(spikes_length)*4,nrow = length(bam_input))

  for (i in 1:length(bam_input)){

    matrix_allele_counts[i,]=rlist::list.cbind(list_all_complete[[i]])
  }

  colnames(matrix_allele_counts)=colnames(rlist::list.cbind(list_all_complete[[1]]))
  row.names(matrix_allele_counts)=cell_names


  name_position_allele=colnames(matrix_allele_counts)
  name_position=rep(0,length(colnames(matrix_allele_counts)))
  for (i in 1:length(colnames(matrix_allele_counts))){
    part_1=strsplit(colnames(matrix_allele_counts)[i],split="_")[[1]][1]
    part_2=strsplit(colnames(matrix_allele_counts)[i],split="_")[[1]][4]
    name_position[i]=paste(part_1,part_2,sep="_")
    print(i)
  }


  name_position_unique=unique(name_position)



  return(list(matrix_allele_counts,colnames(matrix_allele_counts),name_position))

}




#output_SNP_sc_mt_new=Get_raw_counts_allele(bam_input,path_fasta,cell_names)

#output_SNP_sc_mt_new=Get_raw_counts_allele(bam_input,path_fasta,cell_names)

#result=Get_raw_counts_allele(opt$'bam_input',opt$'path_fasta',get(opt$'cell_names'))

start_time <- Sys.time()
result=get_raw_counts_allele(get(opt$'bam_input'),opt$'path_fasta',get(opt$'cell_names'),opt$'cores_number')
end_time <- Sys.time()

print(end_time - start_time)





save(result,file=opt$output)



#setwd("/home/ies/gabriele.lubatti/revision_heteroplasmy/output_R")
#save(output_SNP_sc_mt_new,file="output_SNP_sc_mt_new.Rda")

