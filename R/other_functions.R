get_distribution=function(matrix,FUNCTION,index=NULL){
  base=colnames(matrix)
  distribution=rep(0,length(base))
  if(is.null(index)){
    for (i in 1:length(base)){



      Y=as.vector(matrix[,base[i]])
      if (FUNCTION=='mean'){
        Y=mean(Y)}
      if (FUNCTION=='max'){
        Y=max(Y)}
      if (FUNCTION=='min'){
        Y=min(Y)}
      if (FUNCTION=='median'){
        Y=median(Y)}
      if (FUNCTION=='sum'){
        Y=sum(Y)}

      distribution[i]=Y
      names(distribution)[i]=base[i]
    }}



  else{

    for (i in 1:length(base)){
      index_cell=index[[which(names(index)==base[i])]]


      Y=as.vector(matrix[index_cell,base[i]])
      if (FUNCTION=='mean'){
        Y=mean(Y)}
      if (FUNCTION=='max'){
        Y=max(Y)}
      if (FUNCTION=='min'){
        Y=min(Y)}
      if (FUNCTION=='median'){
        Y=median(Y)}
      if (FUNCTION=='sum'){
        Y=sum(Y)}

      distribution[i]=Y
      names(distribution)[i]=base[i]
    }}
  return(distribution)
}




filter_bases=function(heteroplasmy_matrix,min_heteroplasmy,min_cells,index=NULL){
  position=colnames(heteroplasmy_matrix)
  filter=rep(0,length(position))
  if(is.null(index)){
    for (i in 1:length(position)){
      Y=data.frame(t(heteroplasmy_matrix[,position[i]]))
      filter[i]=sum(Y>min_heteroplasmy)
    }
    heteroplasmy_matrix_select=heteroplasmy_matrix[,which(filter>min_cells)]


  }
  else{
    for (i in 1:length(position)){
      Y=data.frame(t(heteroplasmy_matrix[as.numeric(index[[position[i]]]),position[i]]))
      filter[i]=sum(Y>min_heteroplasmy)
    }
    heteroplasmy_matrix_select=heteroplasmy_matrix[,which(filter>min_cells)]}
  return(colnames(heteroplasmy_matrix_select))
}




detect_insertion=function(ref_sequence,different_sequence,length_comparison=10){
  max_length_insertion=length(different_sequence)-length(ref_sequence)
  for (i in 1:(length(ref_sequence)-length_comparison)){
    if (length(different_sequence)>length(ref_sequence)){
      if(different_sequence[i]!=ref_sequence[i]){
        for(k in 1:max_length_insertion){
          if(all(different_sequence[(i+k):(i+k+length_comparison)]==ref_sequence[i:(i+length_comparison)])){
            # print(k)
            different_sequence=different_sequence[-seq(i,(i+k-1))]
            break
          }
        }
      }}
    else{}
  }
  return(different_sequence)}


















