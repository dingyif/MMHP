clusterAnalysis <- function(start, end, state="social"){
  if(state == "social"){
    if(start == 0){
      input_matrix <- adjm[round(end),,]-utility_state_cumcount[round(end),,]
    }else{
      input_matrix <- adjm[round(end),,]-utility_state_cumcount[round(end),,]-adjm[round(start),,]+utility_state_cumcount[round(start),,]
    }
  }else if(state == "utility"){
    if(start == 0){
      input_matrix <- utility_state_cumcount[round(end),,]
    }else{
      input_matrix <- utility_state_cumcount[round(end),,]-utility_state_cumcount[round(start),,]
    }
  }else if(state == "total"){
    if(start == 0){
      input_matrix <- adjm[round(end),,]
    }else{
      input_matrix <- adjm[round(end),,]-adjm[round(start),,]
    }
  }else{
    return("ERROR")
  }
  
  if(sum(input_matrix)==0){
    return(0)
  }
  
  g_community <- graph_from_adjacency_matrix(input_matrix,mode="undirected")
  cluster_temp <- cluster_leading_eigen(g_community, options=list(maxiter=10000))
  membership <- cluster_temp$membership #cluster_leading_eigen
  number_cluster <- length(unique(membership))
  cluster_member_count <- numeric(number_cluster)
  # cluster_modularity <- modularity(cluster_temp)
  
  within_count <- 0
  between_count <- 0
  for(cluster in c(1:number_cluster)){
    within_count <- within_count + sum(input_matrix[membership==cluster,membership==cluster])
    cluster_member_count[cluster] <- sum(membership==cluster)
    for(another_cluster in c(1:number_cluster)[-cluster]){
      between_count <- between_count + sum(input_matrix[membership==cluster,membership==another_cluster])
    }
  }
  
  between_count_normalizor <- 0
  if(number_cluster == 1){
    between_count_normalizor <- 0
  }else if(number_cluster == 2){
    between_count_normalizor <- cluster_member_count[1]*cluster_member_count[2]
  }else{
    for(cluster in c(1:(number_cluster-1))){
      for(cluster2 in c((cluster+1):number_cluster)){
        between_count_normalizor <- between_count_normalizor + cluster_member_count[cluster]*cluster_member_count[cluster2]
      }
    }
  }
  
  normalized_between_count <- between_count/between_count_normalizor
  normalized_within_count <- within_count/ifelse(sum(sapply(cluster_member_count,choose,k=2))==0,1,sum(sapply(cluster_member_count,choose,k=2)))
  normalized_all_count <- (within_count+between_count)/(between_count_normalizor+sum(sapply(cluster_member_count,choose,k=2)))
  return(normalized_between_count/normalized_within_count)
  ## return(list(normalized_within_count=normalized_within_count,normalized_all_count=normalized_all_count))
  #return(cluster_modularity)
}
